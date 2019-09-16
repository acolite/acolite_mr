## def bundle_test
## lists files in given directory and returns dict with band and file path
## written by Quinten Vanhellemont, RBINS
## 2018-03-12
## modifications: 2018-03-14 (QV) added option to give .tif or metadata.xml files
##                2018-03-14 (QV) improved filtering to include clipped files
##                2018-03-19 (QV) added MS files in filtering

def bundle_test(bundle):
    import os

    if os.path.isdir(bundle):
        files = os.listdir(bundle)
        datafiles = {}
        for i, fname in enumerate(files):
            fn,ext = os.path.splitext(fname)
            if ext == '.json': continue
            if ext not in ['.tif', '.xml']: continue
            band,clp=None,''
            if 'clip' in fn:
                clp='_clip'
                print(fname)
                if 'Analytic_metadata{}.xml'.format(clp) in fname:
                    band = 'metadata'
                if 'AnalyticMS_metadata{}.xml'.format(clp) in fname:
                    band = 'metadata'

                if 'Analytic{}.tif'.format(clp) in fname:
                        band = 'analytic'
                if 'AnalyticMS{}.tif'.format(clp) in fname:
                        band = 'analytic'

                if 'DN_udm{}.tif'.format(clp) in fname:
                    band = 'udm'
                if 'Analytic_SR{}.tif'.format(clp) in fname:
                    band = 'sr'
            else:
                if 'Analytic_metadata.xml' in fname:
                    band = 'metadata'
                if 'AnalyticMS_metadata.xml' in fname:
                    band = 'metadata'
                if 'Analytic.tif' in fname:
                        band = 'analytic'
                if 'AnalyticMS.tif' in fname:
                        band = 'analytic'
                if 'DN_udm.tif' in fname:
                    band = 'udm'
                if 'AnalyticMS_SR.tif' in fname:
                    band = 'sr'

            #else:
            #    band = fn[24:]
            if band is None: continue
            file = '{}/{}'.format(bundle,fname)
            datafiles[band] = {"path":file, "fname":fname}
    else:
        datafiles = {}
        if '.tif' in bundle:
            fname = os.path.basename(bundle)
            dname = os.path.dirname(bundle)
            fn,ext = os.path.splitext(fname)
            #band = fn[24:]
            if 'Analytic' in fname:
                    band = 'analytic'

            file = '{}/{}'.format(dname,fname)
            datafiles[band] = {"path":file, "fname":fname}
            mname = '{}_metadata.xml'.format(fn)
            metafile = '{}/{}'.format(dname,mname)
            if os.path.isfile(metafile):
                datafiles['metadata'] = {"path":metafile, "fname":mname}

        if '_metadata.xml' in bundle:
            fname = os.path.basename(bundle)
            dname = os.path.dirname(bundle)
            fn,ext = os.path.splitext(fname)
            band = 'metadata'
            file = '{}/{}'.format(dname,fname)
            datafiles[band] = {"path":file, "fname":fname}
            tb = '_'.join(fn.split('_')[0:-1])
            tname = '{}.tif'.format(tb)
            tiffile = '{}/{}'.format(dname,tname)
            if os.path.isfile(tiffile):
                datafiles['analytic'] = {"path":tiffile, "fname":tname}

    return(datafiles)
