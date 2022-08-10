def get_G_frame_stuff(gcd_file):
    from icecube import dataio, dataclasses, icetray
    geo = None
    file = dataio.I3File(str(gcd_file))
    while(file.more()):
        frame = file.pop_frame()
        if frame.Stop == icetray.I3Frame.Geometry:
            geo = frame['I3Geometry'].omgeo
        else:
            continue
    if geo is None:
        raise LookupError('Cannot find a G frame in the provided file {}'.format(gcd_file))

    om_x = []
    om_z = []
    strings_x = []
    strings_y = []
    finished = []
    for omkey, omgeo in geo:
        if omgeo.position.z > 1000:
            continue
        om_x.append(omgeo.position.x)
        om_z.append(omgeo.position.z)
        if omkey.om == 1:
            if omkey.string not in finished:
                finished.append(omkey.string)
                strings_x.append(omgeo.position.x)
                strings_y.append(omgeo.position.y)
    om_x = np.asarray(om_x)
    om_z = np.asarray(om_z)
    strings_x = np.asarray(strings_x)
    strings_y = np.asarray(strings_y)
    z_surface = dataclasses.I3Constants.SurfaceElev - dataclasses.I3Constants.OriginElev
    # om_z -= z_surface
    return om_x, om_z, strings_x, strings_y, z_surface

import numpy as np
om_x, om_z, strings_x, strings_y, surface = get_G_frame_stuff("/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2020.Run134142.Pass2_V0.i3.gz")
np.savez('dom_loc.npz', om_x=om_x, om_z=om_z, strings_x = strings_x, strings_y = strings_y, surface=surface)

