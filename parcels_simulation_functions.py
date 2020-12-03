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


def vertical_diffusion_run(k_z, w_10, w_rise, diffusion_type, boundary='Mixed'):
    # Create the fieldset
    fieldset = create_fieldset(k_z=k_z, w_10=w_10, w_rise=w_rise, diffusion_type=diffusion_type)
    # Create the particle set
    pset = ParticleSet(fieldset=fieldset, pclass=JITParticle, lon=[0.5]*SET.p_number, lat=[0.5]*SET.p_number,
                       depth=[0]*SET.p_number)
    output_file = pset.ParticleFile(name=utils.get_parcels_output_name(k_z, w_10, w_rise, diffusion_type, boundary),
                                    outputdt=SET.dt_out)
    # # Determine the particle behavior
    # if diffusion_type == 'Rouse':
    #     kernel = pset.Kernel(simple_vertical_diffusion_Rouse)
    # elif diffusion_type == 'Kukulka':
    #     kernel = pset.Kernel(simple_vertical_diffusion_Kukulka)
    # elif diffusion_type == 'KPP':
    #     kernel = pset.Kernel(simple_vertical_diffusion_KPP)
    # # The actual integration
    # pset.execute(kernel, runtime=SET.runtime, dt=SET.dt_int, output_file=output_file)
    # output_file.export()



def create_fieldset(k_z, w_10, w_rise, diffusion_type):
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
    return fieldset


def simple_vertical_diffusion_mixed_layer_boundary(particle, fieldset, time):
    # According to Ross & Sharples (2004), first the deterministic part of equation 1
    dK_z_p = fieldset.dK_z[time, particle]
    deterministic = dK_z_p * particle.dt

    # The random walk component
    R = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    bz = math.sqrt(fieldset.K_z(time, particle.depth + 0.5 * dK_z_p * particle.dt, particle.lat, particle.lon))

    # Rise velocity component
    rise = fieldset.wrise * particle.dt

    # Total movement
    w_total = deterministic + R * bz + rise

    # Dealing with boundary condition by having a thin mixed layer, and if the particle is not within the mixed layer,
    # then we have simple diffusion according to the equation
    boundary_depth = 0.22  # m
    potential = particle.depth + w_total
    if potential < boundary_depth:
        particle.depth = ParcelsRandom.uniform(0, 1.) * boundary_depth
    else:
        particle.depth = potential

def simple_vertical_diffusion_Rouse(particle, fieldset, time):
    rho_w_p = fieldset.rho_w # density sea water (kg/m^3)
    rho_a_p = fieldset.rho_a # density air (kg/m^3)
    w_10_p = fieldset.w10 # 10 meter wind speed
    C_D = min(max(1.2E-3, 1.0E-3*(0.49+0.065*w_10_p)), 2.12E-3)
    u_s = math.sqrt(rho_a_p / rho_w_p * C_D) * math.fabs(w_10_p) # shear velocity
    k = fieldset.vk # von Karman constant
    K_s = k * u_s * particle.depth + fieldset.Kz
    dK_s = k * u_s

    # According to Ross & Sharples (2004), first the deterministic part of equation 1
    deterministic = 0#dK_s * particle.dt

    # The random walk component
    R = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    bz = math.sqrt(2 * k * u_s * (particle.depth + 0.5 * dK_s * particle.dt) + fieldset.Kz)

    # Rise velocity component
    rise = fieldset.wrise * particle.dt

    # Total movement
    w_total = deterministic + R * bz * rise

    # Dealing with boundary condition by having a thin mixed layer, and if the particle is not within the mixed layer,
    # then we have simple diffusion according to the equation
    boundary_depth = 0.22 # m
    potential = particle.depth + w_total
    if potential < boundary_depth:
        particle.depth = ParcelsRandom.uniform(0, 1.) * boundary_depth
    else:
        particle.depth = potential


def simple_vertical_diffusion_Kukulka(particle, fieldset, time):
    rho_w_p = fieldset.rho_w # density sea water (kg/m^3)
    rho_a_p = fieldset.rho_a # density air (kg/m^3)
    w_10_p = fieldset.w10 # 10 meter wind speed
    wave_age_p = fieldset.wave_age # wave age of developed wave field (Kukulka et al., 2012)
    g_p = fieldset.g # Gravitational acceleration (
    C_D = min(max(1.2E-3, 1.0E-3*(0.49+0.065*w_10_p)), 2.12E-3)
    u_s = math.sqrt(rho_a_p / rho_w_p * C_D) * math.fabs(w_10_p) # shear velocity
    u_air = math.sqrt(C_D) * math.fabs(w_10_p)
    H_s = 0.96 * g_p ** (-1) * wave_age_p ** 1.5 * u_air ** 2  # significant wave height (m)
    k = fieldset.vk # von Karman constant
    K_s = 1.5 * k * u_s * H_s
    dK_s = 0

    # According to Ross & Sharples (2004), first the deterministic part of equation 1
    deterministic = dK_s * particle.dt

    # The random walk component
    R = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    bz = math.sqrt(1.5 * k * u_s * H_s)

    # Rise velocity component
    rise = fieldset.wrise * particle.dt

    # Total movement
    w_total = deterministic + R * bz * rise

    # Dealing with boundary condition by having a thin mixed layer, and if the particle is not within the mixed layer,
    # then we have simple diffusion according to the equation
    boundary_depth = 0.22  # m
    potential = particle.depth + w_total
    if potential < boundary_depth:
        particle.depth = ParcelsRandom.uniform(0, 1.) * boundary_depth
    else:
        particle.depth = potential
    # Check that the particle doesn't go through the surface, and doesn't go below 200 m for convenience
    # if particle.depth + w_total >= 0 and particle.depth + w_total <= 200:
    #     particle.depth += w_total
    # elif particle.depth + w_total < 0:
    #     particle.depth = 0

def simple_vertical_diffusion_KPP(particle, fieldset, time):
    rho_w_p = fieldset.rho_w # density sea water (kg/m^3)
    rho_a_p = fieldset.rho_a # density air (kg/m^3)
    w_10_p = fieldset.w10 # 10 meter wind speed
    wave_age_p = fieldset.wave_age # wave age of developed wave field (Kukulka et al., 2012)
    g_p = fieldset.g # Gravitational acceleration (
    C_D = min(max(1.2E-3, 1.0E-3*(0.49+0.065*w_10_p)), 2.12E-3)
    u_s = math.sqrt(rho_a_p / rho_w_p * C_D) * math.fabs(w_10_p) # shear velocity
    u_air = math.sqrt(C_D) * math.fabs(w_10_p)
    H_s = 0.96 * g_p ** (-1) * wave_age_p ** 1.5 * u_air ** 2  # significant wave height (m)
    k = fieldset.vk # von Karman constant
    MLD_p = fieldset.MLD
    phi_p = fieldset.phi
    z0_p = fieldset.z0

    # Getting the diffusivity profile according to Boufadel et al. (2020)
    # DOI: https://doi.org/10.1029/2019JC015727
    alpha = (k * u_s)/phi_p
    dK_s = alpha/math.pow(MLD_p, 2) * (MLD_p - particle.depth)*(MLD_p - 3 * particle.depth - 2 * z0_p)

    # According to Ross & Sharples (2004), first the deterministic part of equation 1
    deterministic = dK_s * particle.dt

    # The random walk component
    R = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    K_s = alpha * ((particle.depth + 0.5 * dK_s * particle.dt) + z0_p) * (1 - (particle.depth + 0.5 * dK_s * particle.dt)/MLD_p)**2
    bz = math.sqrt(K_s)

    # Rise velocity component
    rise = fieldset.wrise * particle.dt

    # Total movement
    w_total = deterministic + R * bz + rise

    # Dealing with boundary condition by having a thin mixed layer, and if the particle is not within the mixed layer,
    # then we have simple diffusion according to the equation
    boundary_depth = 0.22  # m
    potential = particle.depth + w_total
    if potential < boundary_depth:
        particle.depth = ParcelsRandom.uniform(0, 1.) * boundary_depth
    else:
        particle.depth = potential
    # Check that the particle doesn't go through the surface, and doesn't go below 200 m for convenience
    # if particle.depth + w_total >= 0 and particle.depth + w_total <= 200:
    #     particle.depth += w_total
    # elif particle.depth + w_total < 0:
    #     particle.depth = 0