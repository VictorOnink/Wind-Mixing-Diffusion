"""
We are going to be doing highly idealised Parcels simulations, with the end goal being to develop a parametrization
for wind mixing for use with particle simulations.

"""

from parcels import FieldSet, ParticleSet, JITParticle, Field
from parcels import ParcelsRandom
import math
import numpy as np
import settings as SET
import utils
import scipy.optimize

def vertical_diffusion_run(k_z, w_10, w_rise, diffusion_type, boundary='Mixed'):
    # Create the fieldset
    fieldset = create_fieldset(k_z=k_z, w_10=w_10, w_rise=w_rise, diffusion_type=diffusion_type, boundary=boundary)
    # Create the particle set
    determine_particle_size(w_rise)
    pset = ParticleSet(fieldset=fieldset, pclass=JITParticle, lon=[0.5]*SET.p_number, lat=[0.5]*SET.p_number,
                       depth=[0]*SET.p_number)
    output_file = pset.ParticleFile(name=utils.get_parcels_output_name(k_z, w_10, w_rise, diffusion_type, boundary),
                                    outputdt=SET.dt_out)
    # Determine the particle behavior
    if boundary is 'Mixed':
        kernel = pset.Kernel(simple_vertical_diffusion_mixed_layer_boundary)
    if boundary is 'Zero_Ceiling':
        kernel = pset.Kernel(simple_vertical_diffusion_zero_ceiling_boundary)
    if boundary is 'Reflect':
        kernel = pset.Kernel(simple_vertical_diffusion_reflective_boundary)
    # The actual integration
    pset.execute(kernel, runtime=SET.runtime, dt=SET.dt_int, output_file=output_file)
    output_file.export()


def determine_mixed_layer(k_z, w_10, w_rise, diffusion_type='KPP'):
    def to_optimize(z_t):
        dK = utils.get_vertical_diffusion_gradient_profile(w_10, z_t, diffusion_type)
        dt = SET.dt_int.seconds
        K = utils.get_vertical_diffusion_profile(w_10, k_z, z_t + 0.5 * dK * dt, diffusion_type)
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


def create_fieldset(k_z, w_10, w_rise, diffusion_type, boundary):
    # Creating the lon and lat grids
    lon = np.linspace(0, 1, num=2)
    lat = np.linspace(0, 1, num=2)
    depth = np.linspace(0, 200, num=400)
    time = np.array([0])
    Time, Depth, Lon, Lat = np.meshgrid(time, depth, lon, lat)
    # Creating the velocity fields and adding those to the field sets
    UV = np.zeros(Lon.shape, dtype=np.float32).reshape(len(time), len(depth), len(lat), len(lon))
    Depth = Depth.reshape(len(time), len(depth), len(lat), len(lon))
    U = Field('U', data=UV, depth=depth, lon=lon, lat=lat)
    V = Field('V', data=UV, depth=depth, lon=lon, lat=lat)
    # Getting the diffusivity fields
    K_z = Field('K_z', data=utils.get_vertical_diffusion_profile(w_10, k_z, Depth, diffusion_type),
                depth=depth, lon=lon, lat=lat)
    dK_z = Field('dK_z', data=utils.get_vertical_diffusion_gradient_profile(w_10, Depth, diffusion_type),
                depth=depth, lon=lon, lat=lat)
    fieldset = FieldSet(U, V)
    fieldset.add_field(K_z)
    fieldset.add_field(dK_z)

    # Add a constant for the vertical diffusivity coefficient
    fieldset.add_constant(name='wrise', value=w_rise)
    # Define the random mixing depth for the mixing boundary condition
    if boundary == 'Mixed':
        fieldset.add_constant(name='mix_depth', value=determine_mixed_layer(k_z, w_10, w_rise, diffusion_type))
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
