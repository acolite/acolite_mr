## QV 2019-02-25 runs Pléiades processing
##
## last modification QV 2019-09-16 added Planet data processing
##                   QV 2020-02-25 added ignore_sr_image keyword
##                   QV 2020-05-19 added pressure and elevation arguments
##                   QV 2020-09-15 added export_geotiff argument, currently only for Planet data

def run_acolite_mr():
    ## import sys to parse arguments
    import sys

    ## for evaluation of bool strings
    import distutils.core

    ## fix matplotlib backend to Agg
    ## skip import if --nogfx is given
    #if ('--nogfx' in sys.argv):
    if True:
        import matplotlib
        matplotlib.use("Agg")

    ## import acolite source
    import acolite as ac

    ## ignore numpy errors
    import numpy as np
    olderr = np.seterr(all='ignore')

    import argparse
    parser = argparse.ArgumentParser(description='ACOLITE MR CLI')
    parser.add_argument('--input', help='Main input bundle containing MS data and metadata')
    parser.add_argument('--output', help='Output directory')
    parser.add_argument('--ancillary_data', help='Get ancillary data (default=False)', default=False)

    parser.add_argument('--uoz_default', help='Default ozone value (default=0.3)', default=0.3)
    parser.add_argument('--uwv_default', help='Default water vapour value (default=1.5)', default=1.5)

    parser.add_argument('--sky_correction', help='Do simple sky glint correction (default=True)', default=True)
    parser.add_argument('--dem_pressure', help='Use DEM to derive pressure (default=False)', default=False)

    parser.add_argument('--aer_models', help='Use these models in the selection process, comma separated (1 Continental 2 Maritime 3 Urban default=1,2)', default=(1,2))

    parser.add_argument('--force_band', help='Force fitting to use this band (Blue,Green,Red,NIR; default=None)', default=None)
    parser.add_argument('--fixed_aot', help='Use this AOT (550nm) with the fixed model (default=None)', default=None)
    parser.add_argument('--fixed_model', help='Use this model with the given (1 Continental 2 Maritime 3 Urban, default=2)', default=None)

    parser.add_argument('--dark_spectrum_full_scene', help='Get dark spectrum from full scene (default=True)', default=True)
    parser.add_argument('--limit', help='Limits for processing, comma separated decimal degrees S,W,N,E (default=None)', default=None)

    parser.add_argument('--output_rgb', help='Output PNG RGB composite (default=True)', default=True)
    parser.add_argument('--pan_sharpen_rgb', help='Pan sharpen the RGB composite (default=False)', default=False)

    parser.add_argument('--ignore_sr_image', help='Skip the SR file provided by Planet (default=True)', default=True)

    parser.add_argument('--elevation', help='Scene elevation in meter (default=None)', default=None)
    parser.add_argument('--pressure', help='Scene pressure in hPa (default=None)', default=None)

    parser.add_argument('--export_geotiff', help='Export NetCDF data also as GeoTIFF - currently Planet only (default=False)', default=False)

    args, unknown = parser.parse_known_args()

    if args.input is None:
        print('No input file given.')
        return(1)

    if args.output is None:
        print('No output directory given.')
        return(1)

    if args.limit is not None:
        limit = [float(s) for s in args.limit.split(',')]
        if len(limit) != 4: args.limit=None
        else: args.limit = limit

    if type(args.ancillary_data) == str: args.ancillary_data = bool(distutils.util.strtobool(args.ancillary_data))

    if type(args.sky_correction) == str: args.sky_correction = bool(distutils.util.strtobool(args.sky_correction))
    if type(args.dem_pressure) == str: args.dem_pressure = bool(distutils.util.strtobool(args.dem_pressure))
    if type(args.dark_spectrum_full_scene) == str: args.dark_spectrum_full_scene = bool(distutils.util.strtobool(args.dark_spectrum_full_scene))

    if type(args.output_rgb) == str: args.output_rgb = bool(distutils.util.strtobool(args.output_rgb))
    if type(args.pan_sharpen_rgb) == str: args.pan_sharpen_rgb = bool(distutils.util.strtobool(args.pan_sharpen_rgb))
    if type(args.export_geotiff) == str: args.export_geotiff = bool(distutils.util.strtobool(args.export_geotiff))

    if args.fixed_aot is not None: args.fixed_aot = float(args.fixed_aot)
    args.uoz_default = float(args.uoz_default)
    args.uwv_default = float(args.uwv_default)

    ## use MOD1 and MOD2 (continental and maritime models)
    luts_all=['PONDER-LUT-201704-MOD1-1013mb', 'PONDER-LUT-201704-MOD2-1013mb', 'PONDER-LUT-201704-MOD3-1013mb']
    if type(args.aer_models) is str:
        args.aer_models = [int(a) for a in args.aer_models.split(',')]
    luts = [luts_all[a-1] for a in args.aer_models]

    fixed_lut = luts_all[1]
    if type(args.fixed_model) == str:
         fixed_lut = luts_all[min(len(luts_all)-1,int(args.fixed_model)-1)]

    ## run the processing
    ac.acolite.acolite_mr_ac(args.input, output=args.output, limit=args.limit,
                            ancillary_data=args.ancillary_data,
                            luts=luts, ignore_sr_image=args.ignore_sr_image,
                            uoz=args.uoz_default,
                            uwv=args.uwv_default,
                            map_rgb=args.output_rgb, map_rgb_rhos=args.output_rgb,
                            export_geotiff=args.export_geotiff,
                            pan_sharpen_rgb=args.pan_sharpen_rgb,
                            sky_correction=args.sky_correction,
                            dem_pressure=args.dem_pressure,
                            pressure = args.pressure,
                            elevation = args.elevation,
                            force_band=args.force_band,
                            fixed_aot550=args.fixed_aot, fixed_lut=fixed_lut,
                            dark_spectrum_full_scene=args.dark_spectrum_full_scene)

if __name__ == '__main__':
    print('Launching ACOLITE MR processing.')
    run_acolite_mr()
