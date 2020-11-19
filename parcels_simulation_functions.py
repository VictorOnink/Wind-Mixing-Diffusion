"""
We are going to be doing highly idealised Parcels simulations, with the end goal being to develop a parametrization
for wind mixing for use with particle simulations.

"""

from parcels import FieldSet, ParticleSet, JITParticle, Field
from parcels import ParcelsRandom as random
import math
import numpy as np
import settings as SET


def vertical_diffusion_run(k_z):
    # Create the fieldset
    fieldset = create_fieldset(k_z=k_z)
    # Create the particle set
    pset = ParticleSet(fieldset=fieldset, pclass=JITParticle, lon=[0.5]*100, lat=[0.5]*100, depth=[0]*100)
    filename = SET.output_dir + 'k_z_{}.py'.format(k_z)
    output_file = pset.ParticleFile(name=filename, outputdt=SET.dt_out)
    # Determine the particle behavior
    kernel = pset.Kernel(simple_vertical_diffusion)
    # The actual integration
    pset.execute(kernel, runtime=SET.runtime, dt=SET.dt_int, output_file=output_file)
    output_file.export()


def create_fieldset(k_z):
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
    return fieldset


def simple_vertical_diffusion(particle, fieldset, time):
    dWz = random.uniform(-1., 1.) * math.sqrt(math.fabs(particle.dt) * 3)
    bz = math.sqrt(2 * fieldset.Kz[time, particle.depth, particle.lat, particle.lon])
    # Check that the particle doesn't go through the surface
    if particle.lon + bz * dWz >= 0:
        particle.lon += bz * dWz