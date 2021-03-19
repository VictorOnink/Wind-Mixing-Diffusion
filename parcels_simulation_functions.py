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
import sys


def vertical_diffusion_run(w_10, w_rise, diffusion_type, alpha, boundary='Mixed'):
    """
    General function that runs a parcels simulation for the given parameters
    """
    # Create the fieldset
    fieldset = create_fieldset(w_10=w_10, w_rise=w_rise, diffusion_type=diffusion_type, boundary=boundary, alpha=alpha)

    # The particle size that corresponds to the rise velocity
    utils.determine_particle_size(w_rise)

    # Set the random seed
    ParcelsRandom.seed(settings.seed)
    print(ParcelsRandom.random())
    sys.exit()

    # Determine the type of particle class (dependent on whether we have a Markov 0 or Markov 1 diffusion approach)
    pset = create_pset(fieldset=fieldset, boundary=boundary)

    # Setting the output file
    output_file = pset.ParticleFile(name=utils.get_parcels_output_name(w_10, w_rise, diffusion_type, boundary, alpha),
                                    outputdt=settings.dt_out)

    # Determine the type of boundary condition we use, which in turn relates to which type of diffusion we consider:
    kernel = create_kernel(pset=pset, boundary=boundary)

    # The actual integration
    pset.execute(kernel, runtime=settings.runtime, dt=settings.dt_int, output_file=output_file,
                 recovery={ErrorCode.ErrorOutOfBounds: utils.DeleteParticle})
    output_file.export()


def create_pset(fieldset, boundary):
    xy = 0.5  # position in lon and lat, which is at the center of the domain
    if 'Markov' in boundary:
        pclass = Markov_1_Particle
        # Creating the particle set
        pset = ParticleSet(fieldset=fieldset, pclass=pclass, lon=[xy] * settings.p_number, lat=[xy] * settings.p_number,
                           depth=[settings.p_start_depth] * settings.p_number, w_p=[0.0] * settings.p_number)
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
                          'Reduce_dt': reduce_dt_boundary_condition}
    return pset.Kernel(markov['Markov' in boundary]) + pset.Kernel(boundary_condition[boundary.replace('_Markov', '')])


class Markov_0_Particle(JITParticle):
    """
    Compute the potential new particle depth, to which the boundary conditions are applied in order to
    """
    potential = Variable('potential', initial=0, dtype=np.float32, to_write=False)


class Markov_1_Particle(JITParticle):
    """
    If using Markov-1 diffusion, we need to keep track of the w_total (w' + w_rise) term (Koszalka et al., 2013)
    """
    to_write = False
    w_p = Variable('w_p', initial=attrgetter('w_p'), dtype=np.float32, to_write=to_write)
    potential = Variable('potential', initial=0, dtype=np.float32, to_write=to_write)


def create_fieldset(w_10, w_rise, diffusion_type, boundary, alpha):
    """
    Creating the fieldset that we use for the simulationw
    """
    # Creating the lon and lat grids
    lon = np.linspace(0, 1, num=2)
    lat = np.linspace(0, 1, num=2)
    max_depth = settings.max_depth
    depth = np.linspace(0, max_depth, num=1000)
    time = np.array([0])
    Time, Depth, Lon, Lat = np.meshgrid(time, depth, lon, lat)

    # Creating the velocity fields and adding those to the field sets
    UV = np.zeros(Lon.shape, dtype=np.float32).reshape(len(time), len(depth), len(lat), len(lon))
    Depth = Depth.reshape(len(time), len(depth), len(lat), len(lon))
    U = Field('U', data=UV, depth=depth, lon=lon, lat=lat)
    V = Field('V', data=UV, depth=depth, lon=lon, lat=lat)

    # Getting the diffusivity fields, and the gradient of the diffusivity fields. The gradients are calculated
    # analytically
    K_z = Field('K_z', data=utils.get_vertical_diffusion_profile(w_10, Depth, diffusion_type),
                depth=depth, lon=lon, lat=lat)
    dK_z = Field('dK_z', data=utils.get_vertical_diffusion_gradient_profile(w_10, Depth, diffusion_type),
                 depth=depth, lon=lon, lat=lat)
    fieldset = FieldSet(U, V)
    fieldset.add_field(K_z)
    fieldset.add_field(dK_z)

    # Add a constant for the vertical rise velocity
    fieldset.add_constant(name='wrise', value=w_rise)
    # Add a constant for the maximum depth,
    fieldset.add_constant(name='max_depth', value=max_depth)
    # Define the random mixing depth for the mixing boundary condition
    if 'Mixed' in boundary:
        fieldset.add_constant(name='mix_depth', value=utils.determine_mixed_layer(w_10, w_rise, diffusion_type))
    if 'Reduce_dt' in boundary:
        fieldset.add_constant(name='dt_max', value=settings.dt_int.seconds)
        fieldset.add_constant(name='dt_min', value=settings.dt_int.seconds / 2 ** 2)
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
    Sharples (2004), so before we look at any boundary condition!!!
    """
    # Deterministic transport due to gradient in Kz
    dK = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]
    d_gradient = dK * particle.dt

    # The wiener increment, and then the random component
    R = math.sqrt(math.fabs(particle.dt) * 3)
    dW = ParcelsRandom.uniform(-R, R)
    # Kz = fieldset.K_z[time, particle.depth + 0.5 * dK * particle.dt, particle.lat, particle.lon]
    Kz = fieldset.K_z[time, particle.depth, particle.lat, particle.lon]
    d_random = math.sqrt(2 * Kz) * dW

    # The rise velocity
    d_rise = fieldset.wrise * particle.dt

    # The potential particle position
    particle.potential = particle.depth + (d_gradient + d_random + d_rise)


def markov_1_potential_position(particle, fieldset, time):
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

    # Now, the variance of the turbulent displacements from K_z. For dt < T_l, Kz ~= sig^2 dt, else Kz ~= sig^2 T_l. We
    # add 1e-20 to sig2 to prevent numerical issues when Kz -> 0
    alp = fieldset.alpha
    sig2 = Kz / dt
    dsig2 = dKz / dt

    # The change in particle velocity based on sig2
    d_random = math.sqrt(2 * (1 - alp) * sig2 / dt) * dWz
    d_gradient = dsig2 * dt

    # The particle perturbation velocity at time t + 1
    particle.w_p = alp * particle.w_p + d_random + d_gradient



"""
------------------------------------------------------------------------------------------------------------------------
Applying the various boundary conditions
------------------------------------------------------------------------------------------------------------------------
"""


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
    a = particle.depth


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
    displacement. For the bottom boundary, we apply a reflecting boundary condition
    """
    if particle.potential < fieldset.mix_depth:
        particle.depth = ParcelsRandom.uniform(0, 1.) * fieldset.mix_depth
    if particle.potential > fieldset.max_depth:
        overshoot = particle.potential - fieldset.max_depth
        particle.depth = fieldset.max_depth - overshoot
    else:
        particle.depth = particle.potential
