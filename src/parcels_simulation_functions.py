"""
We are going to be doing highly idealised Parcels simulations, with the end goal being to develop a parametrization
for wind mixing for use with particle simulations.

"""

from parcels import FieldSet, ParticleSet, JITParticle, Field, OperationCode, Variable, ErrorCode
from parcels import ParcelsRandom
from operator import attrgetter
import math
import numpy as np
import settings
import utils


def vertical_diffusion_run(w_10, w_rise, diffusion_type, alpha=0, boundary='Mixed', theta=1.0, wave_roughness=False):
    """
    General function that runs a parcels simulation for the given parameters
    :param w_10: wind speed
    :param w_rise: particle rise velocity
    :param diffusion_type: the type of diffusion profile, either KPP or SWB
    :param alpha: the alpha memory term, which only plays a role if we are running a M-1 simulation. Otherwise it does
                  get used
    :param boundary: Setting the boundary condition. If the boundary string contains "Markov", then we are running a
                     M-1 process
    :param theta: Langmuir circulation amplification factor for KPP mixing
    :param wave_roughness: if True, have surface roughness be wave height dependent
    :return: None
    """
    # Create the fieldset
    fieldset = create_fieldset(w_10=w_10, w_rise=w_rise, diffusion_type=diffusion_type, boundary=boundary, alpha=alpha,
                               theta=theta, wave_roughness=wave_roughness)

    # The particle size that corresponds to the rise velocity, according to Enders et al. (2015)
    utils.determine_particle_size(w_rise)

    # Set the random seed
    ParcelsRandom.seed(settings.seed)

    # Determine the type of particle class (dependent on whether we have a Markov 0 or Markov 1 diffusion approach)
    pset = create_pset(fieldset=fieldset, boundary=boundary)

    # Setting the output file
    output_file = pset.ParticleFile(name=utils.get_parcels_output_name(w_10, w_rise, diffusion_type, boundary, alpha,
                                                                       theta=theta, wave_roughness=wave_roughness),
                                    outputdt=settings.dt_out)

    # Determine the type of boundary condition we use, which in turn relates to which type of diffusion we consider:
    kernel = create_kernel(pset=pset, boundary=boundary)

    # The actual integration, where we first let the model run for 11 hours to reach a steady state. Afterwards, we save
    # the profile every 3 minutes so that can examine the variability over time.
    pset.execute(kernel, runtime=settings.spinup_time, dt=settings.dt_int,
                 recovery={ErrorCode.ErrorOutOfBounds: utils.DeleteParticle})
    pset.execute(kernel, runtime=settings.runtime, dt=settings.dt_int, output_file=output_file,
                 recovery={ErrorCode.ErrorOutOfBounds: utils.DeleteParticle})
    output_file.export()


def create_pset(fieldset, boundary):
    """
    creating the particle set, depending on whether it is a M-0 or M-1 simulation.
    :param fieldset: the fieldset object
    :param boundary: the string setting the boundary condition
    :return: the ParticleSet object
    """
    # position in lon and lat, which is at the center of the domain and which doesn't play a role in the simulation
    xy = 0.5
    # The particleset for a M-1 simulation
    if 'Markov' in boundary:
        pclass = Markov_1_Particle
        # Creating the particle set
        pset = ParticleSet(fieldset=fieldset, pclass=pclass, lon=[xy] * settings.p_number, lat=[xy] * settings.p_number,
                           depth=[settings.p_start_depth] * settings.p_number, w_p=[0.0] * settings.p_number)
    # The particle set for a M-0 simulation
    else:
        pclass = Markov_0_Particle
        # Creating the particle set
        pset = ParticleSet(fieldset=fieldset, pclass=pclass, lon=[xy] * settings.p_number, lat=[xy] * settings.p_number,
                           depth=[settings.p_start_depth] * settings.p_number)
    return pset


def create_kernel(pset, boundary):
    # Creating the vertical transport kernel, selecting between Markov-0 and Markov-1 and adding the selected boundary
    # condition
    markov = {False: markov_0_potential_position, True: markov_1_potential_position}
    boundary_condition = {'Mixed': mixed_layer_boundary_condition, 'Reflect': reflect_boundary_condition,
                          'Reduce_dt': reduce_dt_boundary_condition, 'Ceiling': ceiling_boundary_condition}
    return pset.Kernel(markov['Markov' in boundary]) + pset.Kernel(boundary_condition[boundary.replace('_Markov', '')])


class Markov_0_Particle(JITParticle):
    """
    For a M-0 simulation, we only need to add the potential depth of the particle as variable, which is the depth of the
    particle before we apply the boundary condition
    """
    potential = Variable('potential', initial=0, dtype=np.float32, to_write=False)


class Markov_1_Particle(JITParticle):
    """
    With a M-1 simulation, we also need to keep track of the perturbation velocity in addition to the potential depth.
    """
    to_write = False
    w_p = Variable('w_p', initial=attrgetter('w_p'), dtype=np.float32, to_write=to_write)
    potential = Variable('potential', initial=0, dtype=np.float32, to_write=to_write)


def create_fieldset(w_10, w_rise, diffusion_type, boundary, alpha, theta, wave_roughness):
    """
    Creating the fieldset that we use for the simulationw
    """
    # Creating the lon and lat grids
    lon = np.linspace(0, 1, num=2)
    lat = np.linspace(0, 1, num=2)
    # Setting the depth domain
    max_depth = settings.max_depth
    depth = np.linspace(0, max_depth, num=settings.depth_levels)
    # The initial time is just 0
    time = np.array([0])
    Time, Depth, Lon, Lat = np.meshgrid(time, depth, lon, lat)

    # Creating the UV velocity fields and adding those to the field sets. These are all set at 0 since they don't play
    # a role in the simulation and are only included to prevent error messages from parcels
    UV = np.zeros(Lon.shape, dtype=np.float32).reshape(len(time), len(depth), len(lat), len(lon))
    Depth = Depth.reshape(len(time), len(depth), len(lat), len(lon))
    U = Field('U', data=UV, depth=depth, lon=lon, lat=lat)
    V = Field('V', data=UV, depth=depth, lon=lon, lat=lat)

    # Getting the diffusivity fields, and the gradient of the diffusivity fields. The gradients are calculated
    # analytically
    K_z = Field('K_z', data=utils.get_vertical_diffusion_profile(w_10, Depth, diffusion_type, theta=theta,
                                                                 wave_roughness=wave_roughness),
                depth=depth, lon=lon, lat=lat)
    dK_z = Field('dK_z', data=utils.get_vertical_diffusion_gradient_profile(w_10, Depth, diffusion_type, theta=theta,
                                                                            wave_roughness=wave_roughness),
                 depth=depth, lon=lon, lat=lat)
    fieldset = FieldSet(U, V)
    fieldset.add_field(K_z)
    fieldset.add_field(dK_z)

    # Add a constant for the vertical rise velocity
    fieldset.add_constant(name='wrise', value=w_rise)
    # Add a constant for the maximum depth of the vertical domain
    fieldset.add_constant(name='max_depth', value=max_depth)
    # Define the random mixing depth for the mixing boundary condition
    if 'Mixed' in boundary:
        fieldset.add_constant(name='mix_depth', value=utils.determine_mixed_layer(w_10, w_rise, diffusion_type))
    # Define the max and min integration timestep values for the reduce_dt boundary condition
    if 'Reduce_dt' in boundary:
        fieldset.add_constant(name='dt_max', value=settings.dt_int.seconds)
        fieldset.add_constant(name='dt_min', value=settings.dt_int.seconds / 2 ** 2)
    # If the simulation is a M-1 simulation, defining the alpha value
    if 'Markov' in boundary:
        fieldset.add_constant(name='alpha', value=alpha)

    return fieldset


"""
------------------------------------------------------------------------------------------------------------------------
Calculating the vertical transport due to stochastic transport
------------------------------------------------------------------------------------------------------------------------
"""


def markov_0_potential_position(particle, fieldset, time):
    """
    Here we only calculate the potential particle position following the Markov-0 approach as described by Ross &
    Sharples (2004), so before we look at any boundary condition
    """
    # Deterministic transport due to gradient in Kz
    dK = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]
    d_gradient = dK * particle.dt

    # The wiener increment, and then the random component
    R = math.sqrt(math.fabs(particle.dt) * 3)
    dW = ParcelsRandom.uniform(-R, R)
    Kz = fieldset.K_z[time, particle.depth, particle.lat, particle.lon]
    d_random = math.sqrt(2 * Kz) * dW

    # The rise velocity
    d_rise = fieldset.wrise * particle.dt

    # The potential particle position
    particle.potential = particle.depth + (d_gradient + d_random + d_rise)


def markov_1_potential_position(particle, fieldset, time):
    """
    This calculates the potential particle position following the M-1 approach outlined by equation 3 from Mofakham &
    Ahmadi (2020)
    """
    # Getting the rise (w_r) and perturbation (w_p) velocities
    w_P, w_r = particle.w_p, fieldset.wrise

    # The potential particle depth
    particle.potential = particle.depth + (w_P + w_r) * particle.dt

    # First, the Wiener increment with zero mean and variance = dt
    dt = particle.dt
    R = math.sqrt(math.fabs(particle.dt) * 3)
    dWz = ParcelsRandom.uniform(-R, R)

    # Next, the vertical diffusion term
    Kz = fieldset.K_z[time, particle.depth, particle.lat, particle.lon]
    dKz = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]

    # Now, the variance of the turbulent displacements from K_z. For dt < T_l, Kz ~= sig^2 dt, else Kz ~= sig^2 T_l.
    sig2 = Kz / dt
    dsig2 = dKz / dt

    # The change in particle velocity based on sig2
    d_random = math.sqrt(2 * (1 - fieldset.alpha) * sig2 / dt) * dWz
    d_gradient = dsig2 * dt

    # The particle perturbation velocity at time t + 1
    particle.w_p = fieldset.alpha * particle.w_p + d_random + d_gradient



"""
------------------------------------------------------------------------------------------------------------------------
Applying the various boundary conditions
------------------------------------------------------------------------------------------------------------------------
"""


def ceiling_boundary_condition(particle, fieldset, time):
    """
    Reflecting boundary conditions at the bottom of the model domain (z = fieldset.max_depth), while if a particle
    crosses the surface boundary layer, the particle depth is set at z=0.
    """
    if particle.potential > 0:
        if particle.potential > fieldset.max_depth:
            overshoot = particle.potential - fieldset.max_depth
            particle.depth = fieldset.max_depth - overshoot
        else:
            particle.depth = particle.potential
    else:
        particle.depth = 0


def reflect_boundary_condition(particle, fieldset, time):
    """
    Reflecting boundary conditions at both the ocean surface and at the bottom of the model domain (z = fieldset.max_depth)
    """
    if particle.potential > 0:
        if particle.potential > fieldset.max_depth:
            overshoot = particle.potential - fieldset.max_depth
            particle.depth = fieldset.max_depth - overshoot
        else:
            particle.depth = particle.potential
    else:
        particle.depth = math.fabs(particle.potential)


def reduce_dt_boundary_condition(particle, fieldset, time):
    """
    If the particle crosses the boundary (either surface or bottom), the time step is halved until dt=min_dt. If the
    particle then still tries to cross then boundary, then a reflective boundary condition is applied.

    Note: The approach generally works, but it is very slow
    """
    if particle.potential > 0:
        if particle.potential > fieldset.max_depth:
            particle.dt /= 2
            if particle.dt < fieldset.dt_min:
                overshoot = particle.potential - fieldset.max_depth
                particle.depth = fieldset.max_depth - overshoot
            else:
                return OperationCode.Repeat
        else:
            particle.depth = particle.potential
            if particle.dt < fieldset.dt_max:
                particle.update_next_dt(particle.dt * 2)
    else:
        particle.dt /= 2
        if particle.dt < fieldset.dt_min:
            particle.depth = math.fabs(particle.potential)
        else:
            return OperationCode.Repeat


def mixed_layer_boundary_condition(particle, fieldset, time):
    """
    Following the approach of Ross & Sharples, there is a mixed layer at the surface boundary where it is assumed all
    particles are distributed homogenously, where the depth of this layer is dependent on the maximum possible
    displacement. For the bottom boundary, we apply a reflecting boundary condition.
    """
    if particle.potential < fieldset.mix_depth:
        particle.depth = ParcelsRandom.uniform(0, 1.) * fieldset.mix_depth
    elif particle.potential > fieldset.max_depth:
        overshoot = particle.potential - fieldset.max_depth
        particle.depth = fieldset.max_depth - overshoot
    else:
        particle.depth = particle.potential
