"""
We are going to be doing highly idealised Parcels simulations, with the end goal being to develop a parametrization
for wind mixing for use with particle simulations.

"""

from parcels import FieldSet, ParticleSet, JITParticle, Field, OperationCode, Variable, ErrorCode
from parcels import ParcelsRandom
from operator import attrgetter
import math
import numpy as np
import settings as SET
import utils
import scipy.optimize


def vertical_diffusion_run(w_10, w_rise, diffusion_type, boundary='Mixed'):
    """
    General function that runs a parcels simulation for the given parameters
    """
    # Create the fieldset
    fieldset = create_fieldset(w_10=w_10, w_rise=w_rise, diffusion_type=diffusion_type, boundary=boundary)

    # The particle size that corresponds to the rise velocity
    determine_particle_size(w_rise)

    # Set the random seed
    ParcelsRandom.seed(SET.seed)

    # Determine the type of particle class (dependent on whether we have a Markov 0 or Markov 1 diffusion approach)
    pset = create_pset(fieldset=fieldset, w_rise=w_rise, boundary=boundary)

    # Setting the output file
    output_file = pset.ParticleFile(name=utils.get_parcels_output_name(w_10, w_rise, diffusion_type, boundary),
                                    outputdt=SET.dt_out)

    # Determine the type of boundary condition we use, which in turn relates to which type of diffusion we consider:
    kernel = create_kernel(pset=pset, boundary=boundary)

    # The actual integration
    pset.execute(kernel, runtime=SET.runtime, dt=SET.dt_int, output_file=output_file,
                 recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
    output_file.export()


def create_pset(fieldset, w_rise, boundary):
    xy = 0.5  # position in lon and lat, which is at the center of the domain
    if 'Markov' in boundary:
        pclass = Markov_1_Particle
        initial_w_total = w_rise + (np.random.rand(SET.p_number) - 0.5) * 2 * SET.w_prime
        # Creating the particle set
        pset = ParticleSet(fieldset=fieldset, pclass=pclass, lon=[xy] * SET.p_number, lat=[xy] * SET.p_number,
                           depth=[SET.p_start_depth] * SET.p_number, w_total=initial_w_total)
    else:
        pclass = Markov_0_Particle
        # Creating the particle set
        pset = ParticleSet(fieldset=fieldset, pclass=pclass, lon=[xy] * SET.p_number, lat=[xy] * SET.p_number,
                           depth=[SET.p_start_depth] * SET.p_number)
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
    If using Markov-1 diffusion, we need to keep track of the w_total (w' + w_rise) term (Koszalka et al., 2013)
    """
    potential = Variable('potential', initial=0, dtype=np.float32, to_write=False)


class Markov_1_Particle(JITParticle):
    """
    If using Markov-1 diffusion, we need to keep track of the w_total (w' + w_rise) term (Koszalka et al., 2013)
    """
    w_total = Variable('w_total', initial=attrgetter('w_total'), dtype=np.float32, to_write=False)
    potential = Variable('potential', initial=0, dtype=np.float32, to_write=False)


def lagrangian_integral_timescale(w_10):
    """
    Determining the Lagrangian integral timescale following Denman & Gargett (1983). If the mixed layer depth is deeper
    than the turbulent Ekman layer thickness, then the Lagrangian integral timescale is dependent only on the coriolis
    parameter. Otherwise, when the mixed layer depth is less than the Ekman layer thickness, then the timescale becomes
    inversely proportional to the surface wind speed

    The ekman timescale is approximately 30 minutes at latitude 45 (MLD > L_E)
    """
    # Friction velocity
    u_w = math.sqrt(utils.determine_tau(w_10, SET.rho_a) / SET.rho_w)
    # RMS eddy velocity
    u_t = 2 * u_w
    # Coriolis parameter
    f = 2 * 7.2921e-5 * math.sin(SET.latitude)
    # The turbulent Ekman layer thickness
    L_e = 0.4 * u_w / f
    # Now determining the time scale
    if L_e < SET.MLD:
        T_L = 0.2 / f
    else:
        T_L = SET.MLD / u_t
    print("The lagrangian integral timescale is {} minutes".format(T_L / 60.))
    return SET.dt_int.seconds  # SET.dt_int.seconds * 2 * 15


def determine_mixed_layer(w_10, w_rise, diffusion_type='KPP'):
    """
    Determining the depth of the surface layer within particles are assumed to be distributed homogeneously, following
    the boundary condition approach of Ross & Sharples (2004)
    """

    def to_optimize(z_t):
        dK = utils.get_vertical_diffusion_gradient_profile(w_10, np.array([z_t]), diffusion_type)
        dt = SET.dt_int.seconds
        K = utils.get_vertical_diffusion_profile(w_10, np.array([z_t + 0.5 * dK * dt]), diffusion_type)
        RHS = dK * dt + np.sqrt(6 * K * dt) + w_rise * dt
        return np.abs(z_t - RHS)

    mixing_depth = scipy.optimize.minimize_scalar(to_optimize, bounds=[0, 100], method='bounded').x
    print('The surface turbulent mixed layer depth {} m'.format(mixing_depth))
    return mixing_depth


def determine_particle_size(w_rise):
    """
    Determining the size of an elliptical particle that corresponds to the rise velocity, following the approach
    of Poulain et al. (2019)
    """

    def Re(L):
        # Return the Reynolds number
        return L * np.abs(w_rise) / (SET.mu / SET.rho_w)

    def to_optimize(L):
        Reynold = Re(L)
        left = 240. / (np.pi * Reynold) * (1 + 0.138 * Reynold ** 0.792) * w_rise ** 2
        right = 2. / 15 * L * (SET.rho_p / SET.rho_w - 1) * SET.g
        return np.abs(left - right)

    particle_size = scipy.optimize.minimize_scalar(to_optimize, bounds=[0, 100], method='bounded').x
    print('The particle size according to Poulain et al. (2019) is {} m'.format(particle_size))


def create_fieldset(w_10, w_rise, diffusion_type, boundary):
    """
    Creating the fieldset that we use for the simulationw
    """
    # Creating the lon and lat grids
    lon = np.linspace(0, 1, num=2)
    lat = np.linspace(0, 1, num=2)
    max_depth = SET.max_depth
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
        fieldset.add_constant(name='mix_depth', value=determine_mixed_layer(w_10, w_rise, diffusion_type))
    if 'Reduce_dt' in boundary:
        fieldset.add_constant(name='dt_max', value=SET.dt_int.seconds)
        fieldset.add_constant(name='dt_min', value=SET.dt_int.seconds / 2 ** 2)
    if 'Markov' in boundary:
        fieldset.add_constant(name='T_L', value=lagrangian_integral_timescale(w_10))

    return fieldset


def markov_0_potential_position(particle, fieldset, time):
    """
    Here we only calculate the potential particle position following the Markov-0 approach as described by Ross &
    Sharples (2004), so before we look at any boundary condition!!!
    """
    # Deterministic transport due to gradient in Kz
    dK = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]
    d_gradient = dK * particle.dt

    # The wiener increment, and then the random component
    dW = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    d_random = math.sqrt(2 * fieldset.K_z[time, particle.depth + 0.5 * dK * particle.dt, particle.lat, particle.lon]) * dW

    # The rise velocity
    d_rise = fieldset.wrise * particle.dt

    # The potential particle position
    particle.potential = particle.depth + (d_gradient + d_random + d_rise)


def markov_1_potential_position(particle, fieldset, time):
    """
    Here we only calculate the potential particle position following the Markov-1 approach as described by following
    eq. 4 of Brickman & Smith (2002), so before we look at any boundary condition!!!

    We assume a constant w_rise and stationary turbulence
    """
    # First, the Wiener increment with zero mean and variance = dt
    dt = particle.dt
    dWz = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(dt)))

    # Next, the vertical diffusion term, and the gradient of the vertical diffusion
    Kz = fieldset.K_z[time, particle.depth, particle.lat, particle.lon]
    dKz = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]

    # Now, the variance of the turbulent displacements from K_z. For dt < T_l, Kz ~= sig^2 dt, else Kz ~= sig^2 T_l
    T_l = fieldset.T_L
    sig2 = max(Kz / dt, Kz / T_l)
    dsig2 = max(dKz / dt, dKz / T_l)

    # Getting the total (w_T), rise (w_r) and perturbation (w_p) velocities, where w_T = w_r + w_p
    w_T, w_r = particle.w_total, fieldset.wrise
    w_p = w_T - w_r

    # The change in particle velocity
    d_history = - (w_p / T_l) * dt
    d_gradient = 0.5 * dsig2 * dt * (1 + w_p * (w_p + w_r) / sig2)
    d_random = math.sqrt(2 * sig2 / T_l) * dWz
    particle.w_total += (d_history + d_gradient + d_random)

    # The potential particle depth
    particle.potential = particle.depth + particle.w_total * particle.dt


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
    displacement. For the bottom boundary, we apply a reflecting boundary condition
    """
    if particle.potential < fieldset.mix_depth:
        particle.depth = ParcelsRandom.uniform(0, 1.) * fieldset.mix_depth
    if particle.potential > fieldset.max_depth:
        overshoot = particle.potential - fieldset.max_depth
        particle.depth = fieldset.max_depth - overshoot
    else:
        particle.depth = particle.potential


def DeleteParticle(particle, fieldset, time):
    particle.delete()