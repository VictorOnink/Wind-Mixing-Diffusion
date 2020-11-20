"""
We are going to be doing highly idealised Parcels simulations, with the end goal being to develop a parametrization
for wind mixing for use with particle simulations.

"""

from parcels import FieldSet, ParticleSet, JITParticle, Field
from parcels import ParcelsRandom
import math
import numpy as np
import settings as SET


def vertical_diffusion_run(k_z, w_10, w_rise):
    # Create the fieldset
    fieldset = create_fieldset(k_z=k_z, w_10=w_10, w_rise=w_rise)
    # Create the particle set
    pset = ParticleSet(fieldset=fieldset, pclass=JITParticle, lon=[0.5]*SET.p_number, lat=[0.5]*SET.p_number,
                       depth=[0]*SET.p_number)
    output_file = pset.ParticleFile(name=get_parcels_output_name(k_z, w_10, w_rise), outputdt=SET.dt_out)
    # Determine the particle behavior
    kernel = pset.Kernel(simple_vertical_diffusion)
    # The actual integration
    pset.execute(kernel, runtime=SET.runtime, dt=SET.dt_int, output_file=output_file)
    output_file.export()


def create_fieldset(k_z, w_10, w_rise):
    # Creating the lon and lat grids
    lon = np.linspace(0, 1, num=2)
    lat = np.linspace(0, 1, num=2)
    Lon, Lat = np.meshgrid(lon,lat)
    # Creating the velocity fields and adding those to the field sets
    UV = np.zeros(Lon.shape, dtype=np.float32)
    U = Field('U', data=UV, lon=lon, lat=lat)
    V = Field('V', data=UV, lon=lon, lat=lat)
    fieldset = FieldSet(U, V)
    # Add a constant for the vertical diffusivity coefficient
    fieldset.add_constant(name='Kz', value=k_z)
    fieldset.add_constant(name='w10', value=w_10)
    fieldset.add_constant(name='wrise', value=w_rise)
    return fieldset


def simple_vertical_diffusion(particle, fieldset, time):
    rho_w = 1027 # density sea water (kg/m^3)
    rho_a = 1.22 # density air (kg/m^3)
    w_10 = fieldset.w10 # 10 meter wind speed
    C_D = min(max(1.2E-3, 1.0E-3*(0.49+0.065*w_10)), 2.12E-3)
    u_s = math.sqrt(rho_a / rho_w * C_D) * math.fabs(w_10) # shear velocity
    k = 0.4 # von Karman constant
    K_s = k * u_s * particle.depth + fieldset.Kz

    # Rise velocity
    w_rise = fieldset.wrise
    # The random walk
    dWz = ParcelsRandom.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    bz = math.sqrt(2 * K_s)
    # Total movement
    w_total = bz * dWz + w_rise

    # Check that the particle doesn't go through the surface, and doesn't go below 200 m for convenience
    if particle.depth + w_total >= 0 and particle.depth + w_total <= 200:
        particle.depth += w_total


def get_parcels_output_name(k_z, w_10, w_rise):
    return SET.output_dir + 'k_z_{}_w10_{}_w_rise_{}.nc'.format(k_z, w_10, w_rise)