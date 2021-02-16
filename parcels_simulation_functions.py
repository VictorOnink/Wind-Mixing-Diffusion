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
    if 'Markov' in boundary:
        pclass = MarkovParticle
        initial_w_total = w_rise + (np.random.rand(SET.p_number) - 1) * 2 * SET.w_prime
        # Creating the particle set
        pset = ParticleSet(fieldset=fieldset, pclass=pclass, lon=[0.5] * SET.p_number, lat=[0.5] * SET.p_number,
                           depth=[SET.p_start_depth] * SET.p_number, w_total=initial_w_total)
    else:
        pclass = JITParticle
        # Creating the particle set
        pset = ParticleSet(fieldset=fieldset, pclass=pclass, lon=[0.5] * SET.p_number, lat=[0.5] * SET.p_number,
                           depth=[SET.p_start_depth] * SET.p_number)

    # Setting the output file
    output_file = pset.ParticleFile(name=utils.get_parcels_output_name(w_10, w_rise, diffusion_type, boundary),
                                    outputdt=SET.dt_out)

    # Determine the type of boundary condition we use, which in turn relates to which type of diffusion we consider:
    boundary_dict = {'Mixed': markov_0_mixed_layer_boundary, 'Reflect': markov_0_reflect,
                     'Reduce_dt': markov_0_reduce_dt, 'Mixed_Markov': markov_1_mixed_layer_boundary,
                     'Reflect_Markov': markov_1_reflect, 'Reduce_dt_Markov': markov_1_reduce_dt}
    kernel = pset.Kernel(boundary_dict[boundary])

    # The actual integration
    pset.execute(kernel, runtime=SET.runtime, dt=SET.dt_int, output_file=output_file,
                 recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})
    output_file.export()


class MarkovParticle(JITParticle):
    """
    If using Markov-1 diffusion, we need to keep track of the w_total (w' + w_rise) term (Koszalka et al., 2013)
    """
    w_total = Variable('w_total', initial=attrgetter('w_total'), dtype=np.float32, to_write=False)


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
    return SET.dt_int.seconds * 2 * 15


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


def markov_0_mixed_layer_boundary(particle, fieldset, time):
    """
    Following the approach of Ross & Sharples, there is a mixed layer at the boundary where it is assumed all particles
    are distributed homogenously, where the depth of this layer is dependent on the maximum possible displacement

    Basically, if the new particle position is within the mixed layer (or above it), then the particle is distributed
    randomly throughout this mixing layer (and the particle will never fly through the ocean surface)
    """
    # According to Ross & Sharples (2004), first the deterministic part of equation 1
    dK_z_p = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]
    deterministic = dK_z_p * particle.dt

    # The random walk component
    R = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    bz = math.sqrt(fieldset.K_z[time, particle.depth + 0.5 * dK_z_p * particle.dt, particle.lat, particle.lon])

    # Rise velocity component
    rise = fieldset.wrise * particle.dt

    # Total movement
    w_total = deterministic + R * bz + rise

    # Dealing with boundary condition by having a thin mixed layer, and if the particle is not within the mixed layer,
    # then we have simple diffusion according to the equation
    potential = particle.depth + w_total
    if potential < fieldset.mix_depth:
        particle.depth = ParcelsRandom.uniform(0, 1.) * fieldset.mix_depth
    else:
        particle.depth = potential


def markov_0_reflect(particle, fieldset, time):
    """
    If a particle tries to cross the boundary, then it is reflected back
    """
    # According to Ross & Sharples (2004), first the deterministic part of equation 1
    dK_z_p = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]
    deterministic = dK_z_p * particle.dt

    # The random walk component
    R = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    bz = math.sqrt(fieldset.K_z[time, particle.depth + 0.5 * dK_z_p * particle.dt, particle.lat, particle.lon])

    # Rise velocity component
    rise = fieldset.wrise * particle.dt

    # Total movement
    w_total = deterministic + R * bz + rise

    # The ocean surface acts as a lid off of which the plastic bounces if tries to cross the ocean surface
    potential = particle.depth + w_total
    if potential < 0:
        particle.depth = -1 * potential
    else:
        particle.depth = potential


def markov_0_reduce_dt(particle, fieldset, time):
    """
    If the particle crosses the boundary, the time step is halved until dt=min_dt. If the particle then still tries
    to cross then boundary, then a reflective boundary condition is applied.

    Note: The approach generally works, but it is very slow
    """
    # According to Ross & Sharples (2004), first the deterministic part of equation 1
    dK_z_p = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]
    deterministic = dK_z_p * particle.dt

    # The random walk component
    R = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    bz = math.sqrt(fieldset.K_z[time, particle.depth + 0.5 * dK_z_p * particle.dt, particle.lat, particle.lon])

    # Rise velocity component
    rise = fieldset.wrise * particle.dt

    # Total movement
    w_total = deterministic + R * bz + rise

    # The ocean surface acts as a lid off of which the plastic bounces if tries to cross the ocean surface
    potential = particle.depth + w_total

    if potential > 0:
        particle.depth = potential
        if particle.dt < fieldset.dt_max:
            particle.update_next_dt(particle.dt * 2)
    else:
        particle.dt /= 2
        if particle.dt < fieldset.dt_min:
            particle.depth = -1 * potential
        else:
            return OperationCode.Repeat


def markov_1_mixed_layer_boundary(particle, fieldset, time):
    """
    Following the approach of Ross & Sharples, there is a mixed layer at the boundary where it is assumed all particles
    are distributed homogenously, where the depth of this layer is dependent on the maximum possible displacement

    Basically, if the new particle position is within the mixed layer (or above it), then the particle is distributed
    randomly throughout this mixing layer (and the particle will never fly through the ocean surface)
    """

    # First, the Wiener increment with zero mean and variance = dt
    dt = particle.dt
    dWz = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(dt)))
    # Next, the vertical diffusion term, and the gradient of the vertical diffusion
    Kz = fieldset.K_z[time, particle.depth, particle.lat, particle.lon]
    dKz = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]

    # The new random velocity perturbation, following Koszalka et al. (2013)
    # dw_prime = 1. / fieldset.T_L * (math.sqrt(2 * Kz) * dWz - particle.w_prime * particle.dt)
    # The new random velocity perturbation, following eq. 5 of Brickman & Smith (2002)
    w = particle.w_total
    dw_prime = 1. / fieldset.T_L * (
            (-1 * w + 0.5 * dKz * (w ** 2 * fieldset.T_L / Kz + 1)) * dt + math.sqrt(2 * Kz) * dWz)
    particle.w_total += dw_prime

    # The ocean surface acts as a lid off of which the plastic bounces if tries to cross the ocean surface
    particle.depth += (particle.w_total + fieldset.wrise) * particle.dt

    # Dealing with boundary condition by having a thin mixed layer, and if the particle is not within the mixed layer,
    # then we have simple diffusion according to the equation
    if particle.depth < fieldset.mix_depth:
        particle.depth = ParcelsRandom.uniform(0, 1.) * fieldset.mix_depth
    # If a particle would go deeper than max_depth, it needs to bounce back uo
    if particle.depth > fieldset.max_depth:
        overshoot = particle.depth - fieldset.max_depth
        particle.depth = fieldset.max_depth - overshoot


def markov_1_reflect(particle, fieldset, time):
    # First, the Wiener increment with zero mean and variance = dt
    dt, T_l = particle.dt, fieldset.T_L
    dWz = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(dt)))

    # Next, the vertical diffusion term, and the gradient of the vertical diffusion
    Kz = fieldset.K_z[time, particle.depth, particle.lat, particle.lon]
    dKz = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]

    # Now, the variance of the turbulent displacements from K_z. For dt < T_l, Kz ~= sig^2 dt, else Kz ~= sig^2 T_l
    sig2 = max(Kz / dt, Kz / T_l)
    dsig2 = max(dKz / dt, dKz / T_l)

    # Getting the total (w_T), rise (w_r) and perturbation (w_p) velocities, where w_T = w_r + w_p
    w_T, w_r = particle.w_total, fieldset.wrise
    w_p = w_T - w_r

    # The change in particle velocity, following eq. 4 of Brickman & Smith (2002) but assuming constant w_rise and
    # stationary turbulence. The various terms are:
    d_history = - (w_p / T_l) * dt
    d_gradient = 0.5 * dsig2 * dt * (1 + w_p * (w_p + w_r) / sig2)
    d_random = math.sqrt(2 * sig2 / T_l) * dWz
    particle.w_total += (d_history + d_gradient + d_random)

    # The ocean surface acts as a lid off of which the plastic bounces if tries to cross the ocean surface
    particle.depth = math.fabs(particle.depth + particle.w_total * particle.dt)

    # If a particle would go deeper than max_depth, it needs to bounce back uo
    if particle.depth > fieldset.max_depth:
        overshoot = particle.depth - fieldset.max_depth
        particle.depth = fieldset.max_depth - overshoot


def markov_1_reduce_dt(particle, fieldset, time):
    """
    If the particle crosses the boundary, the time step is halved until dt=min_dt. If the particle then still tries
    to cross then boundary, then a reflective boundary condition is applied.

    Note: The approach generally works, but it is very slow
    """
    # First, the Wiener increment with zero mean and variance = dt
    dt = particle.dt
    dWz = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(dt)))
    # Next, the vertical diffusion term, and the gradient of the vertical diffusion
    Kz = fieldset.K_z[time, particle.depth, particle.lat, particle.lon]
    dKz = fieldset.dK_z[time, particle.depth, particle.lat, particle.lon]

    # The new random velocity perturbation, following Koszalka et al. (2013)
    # dw_prime = 1. / fieldset.T_L * (math.sqrt(2 * Kz) * dWz - particle.w_prime * particle.dt)
    # The new random velocity perturbation, following eq. 5 of Brickman & Smith (2002)
    w = particle.w_total
    dw_total = 1. / fieldset.T_L * (
            (-1 * w + 0.5 * dKz * (w ** 2 * fieldset.T_L / Kz + 1)) * dt + math.sqrt(2 * Kz) * dWz)
    particle.w_total += dw_total
    # The ocean surface acts as a lid off of which the plastic bounces if tries to cross the ocean surface
    potential = particle.depth + (particle.w_total + fieldset.wrise) * particle.dt

    if potential > 0:
        particle.depth = potential
        if particle.dt < fieldset.dt_max:
            particle.update_next_dt(particle.dt * 2)
    else:
        particle.dt /= 2
        if particle.dt < fieldset.dt_min:
            particle.depth = -1 * potential
        else:
            return OperationCode.Repeat


def DeleteParticle(particle, fieldset, time):
    particle.delete()
