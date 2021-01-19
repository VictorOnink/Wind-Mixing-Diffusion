"""
We are going to be doing highly idealised Parcels simulations, with the end goal being to develop a parametrization
for wind mixing for use with particle simulations.

"""

from parcels import FieldSet, ParticleSet, JITParticle, Field, OperationCode, Variable
from parcels import ParcelsRandom
import math
import numpy as np
import settings as SET
import utils
import scipy.optimize
from datetime import timedelta

def vertical_diffusion_run(w_10, w_rise, diffusion_type, boundary='Mixed'):
    # Create the fieldset
    fieldset = create_fieldset(w_10=w_10, w_rise=w_rise, diffusion_type=diffusion_type, boundary=boundary)
    # Create the particle set
    determine_particle_size(w_rise)
    ParcelsRandom.seed(SET.seed)
    if boundary is 'Reflect_Markov':
        pclass = MarkovParticle
    else:
        pclass = JITParticle
    pset = ParticleSet(fieldset=fieldset, pclass=pclass, lon=[0.5]*SET.p_number, lat=[0.5]*SET.p_number,
                       depth=[SET.p_start_depth]*SET.p_number)
    output_file = pset.ParticleFile(name=utils.get_parcels_output_name(w_10, w_rise, diffusion_type, boundary),
                                    outputdt=SET.dt_out)
    # Determine the particle behavior
    if boundary is 'Mixed':
        kernel = pset.Kernel(simple_vertical_diffusion_mixed_layer_boundary)
    elif boundary is 'Zero_Ceiling':
        kernel = pset.Kernel(simple_vertical_diffusion_zero_ceiling_boundary)
    elif boundary is 'Reflect':
        kernel = pset.Kernel(simple_vertical_diffusion_reflective_boundary)
    elif boundary is 'Reduce_dt':
        kernel = pset.Kernel(simple_vertical_diffusion_reduce_dt_boundary)
    elif boundary is 'Reflect_Markov':
        kernel = pset.Kernel(markov1_vertical_diffusion_reflective_boundary)
    # The actual integration
    pset.execute(kernel, runtime=SET.runtime, dt=SET.dt_int, output_file=output_file)
    output_file.export()


class MarkovParticle(JITParticle):
    w_prime = Variable('w_prime', initial=0.0, dtype=np.float32, to_write=False)


def lagrangian_integral_timescale(w_10):
    # Based on Denman & Gargett (1983)
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
    return T_L


def determine_mixed_layer(w_10, w_rise, diffusion_type='KPP'):
    def to_optimize(z_t):
        dK = utils.get_vertical_diffusion_gradient_profile(w_10, z_t, diffusion_type)
        dt = SET.dt_int.seconds
        K = utils.get_vertical_diffusion_profile(w_10, z_t + 0.5 * dK * dt, diffusion_type)
        RHS = dK*dt + np.sqrt(6 * K * dt) + w_rise * dt
        return np.abs(z_t - RHS)
    mixing_depth = scipy.optimize.minimize_scalar(to_optimize, bounds=[0, 100], method='bounded').x
    print('The surface turbulent mixed layer depth {} m'.format(mixing_depth))
    return mixing_depth


def determine_particle_size(w_rise):
    def Re(L):
        return L * np.abs(w_rise) / (SET.mu / SET.rho_w)

    def to_optimize(L):
        Reynold = Re(L)
        left = 240. / (np.pi * Reynold) * (1 + 0.138 * Reynold ** 0.792) * w_rise ** 2
        right = 2. / 15 * L * (SET.rho_p / SET.rho_w - 1) * SET.g
        return np.abs(left - right)

    particle_size = scipy.optimize.minimize_scalar(to_optimize, bounds=[0, 100], method='bounded').x
    print('The particle size according to Poulain et al. (2019) is {} m'.format(particle_size))


def create_fieldset(w_10, w_rise, diffusion_type, boundary):
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
    # Getting the diffusivity fields
    K_z = Field('K_z', data=utils.get_vertical_diffusion_profile(w_10, Depth, diffusion_type),
                depth=depth, lon=lon, lat=lat)
    dK_z = Field('dK_z', data=utils.get_vertical_diffusion_gradient_profile(w_10, Depth, diffusion_type),
                depth=depth, lon=lon, lat=lat)
    fieldset = FieldSet(U, V)
    fieldset.add_field(K_z)
    fieldset.add_field(dK_z)

    # Add a constant for the vertical diffusivity coefficient
    fieldset.add_constant(name='wrise', value=w_rise)
    # Add a constant for the maximum depth, which we set to be the MLD
    fieldset.add_constant(name='max_depth', value=max_depth)
    # Define the random mixing depth for the mixing boundary condition
    if boundary is 'Mixed':
        fieldset.add_constant(name='mix_depth', value=determine_mixed_layer(w_10, w_rise, diffusion_type))
    if boundary is 'Reduce_dt':
        fieldset.add_constant(name='dt_max', value=SET.dt_int.seconds)
        fieldset.add_constant(name='dt_min', value=SET.dt_int.seconds / 2 ** 3)
    if boundary is 'Reflect_Markov':
        fieldset.add_constant(name='T_L', value=lagrangian_integral_timescale(w_10))

    return fieldset


def simple_vertical_diffusion_mixed_layer_boundary(particle, fieldset, time):
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


def simple_vertical_diffusion_zero_ceiling_boundary(particle, fieldset, time):
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

    # If the particle crosses the ocean surface, it is just placed back at the ocean surface
    potential = particle.depth + w_total
    if potential < 0:
        particle.depth = 0
    else:
        particle.depth = potential


def simple_vertical_diffusion_reflective_boundary(particle, fieldset, time):
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


def simple_vertical_diffusion_reduce_dt_boundary(particle, fieldset, time):
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


def markov1_vertical_diffusion_reflective_boundary(particle, fieldset, time):
    # First, the Wiener increment with zero mean and std = dt
    dWz = ParcelsRandom.normalvariate(0, math.sqrt(math.fabs(particle.dt)))
    # Next, the vertical diffusion term, and then the total stochastic term
    Kz = fieldset.K_z[time, particle.depth, particle.lat, particle.lon]

    # The new random velocity perturbation, following Koszalka et al. (2013)
    dw_prime = 1. / fieldset.T_L * (math.sqrt(2 * Kz) * dWz - particle.w_prime * particle.dt)
    particle.w_prime += dw_prime

    # The ocean surface acts as a lid off of which the plastic bounces if tries to cross the ocean surface
    particle.depth = math.fabs(particle.depth + (particle.w_prime + fieldset.wrise) * particle.dt)
    if particle.depth > fieldset.max_depth:
        overshoot = particle.depth - fieldset.max_depth
        particle.depth = fieldset.max_depth - overshoot
