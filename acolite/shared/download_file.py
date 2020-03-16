## def download_file
## download_file with authorisation option
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-14
## modifications: 2018-11-19 (QV) added verbosity option, removed the parallel download
##                2020-04-24 (QV) earthdata login

def download_file(url, file, auth=None, session=None, parallel=False, verbosity=0):

    import requests
    import time

    import os
    file_path = os.path.abspath(file)
    start = time.time()

    if 'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/' in url:
        username = os.environ['EARTHDATA_u']
        password = os.environ['EARTHDATA_p']
        auth = (username, password)

    with requests.Session() as session:
            r1 = session.request('get', url)
            r = session.get(r1.url, auth=auth)

            if (r.ok):
                with open(file_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1024*1024):
                        if chunk: # filter out keep-alive new chunks
                            f.write(chunk)
            else:
                if verbosity > 2: print(r.text)
                raise Exception("File download failed")

    if verbosity > 1:
        print("Downloaded {}, elapsed Time: {:.1f}s".format(url, time.time() - start))

