## def write_rgb
## writes RGB image
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07
## modifications: QV 2017-01-18 converted to function
##                QV 2018-06-14 added mask option

def write_rgb(outputfile, red, green, blue, mask = None, mask_color=(255,255,255),
            rgb_autorange=False, rgb_percentiles=[5,95],
            rgb_range_min=[0.0,0.0,0.0],rgb_range_max=[0.15,0.15,0.15], return_image=False):
    
    from numpy import nanpercentile,dstack
    from acolite import datascl
    from PIL import Image
    
    if rgb_autorange is True:
        rgb_percentiles_r = nanpercentile(red, rgb_percentiles)
        rgb_percentiles_g = nanpercentile(green, rgb_percentiles)
        rgb_percentiles_b = nanpercentile(blue, rgb_percentiles)
        rgb_range_min = [rgb_percentiles_r[0],rgb_percentiles_g[0],rgb_percentiles_b[0]]
        rgb_range_max = [rgb_percentiles_r[1],rgb_percentiles_g[1],rgb_percentiles_b[1]]
        print(rgb_range_min)
        print(rgb_range_max)

    ## rescale reflectances for mapping / RGB composite
    redsc = datascl(red, dmin=rgb_range_min[0], dmax=rgb_range_max[0])
    greensc = datascl(green, dmin=rgb_range_min[1], dmax=rgb_range_max[1])
    bluesc = datascl(blue, dmin=rgb_range_min[2], dmax=rgb_range_max[2])
    #nirsc = datascl(nir, dmin=rgb_range_min[-1], dmax=rgb_range_max[-1])

    if mask is not None:
        redsc[mask != 0] = mask_color[0]
        greensc[mask != 0] = mask_color[1]
        bluesc[mask != 0] = mask_color[2]

    ## stack scaled bands for RGB image
    image = dstack((redsc,greensc,bluesc))

    ## output image    
    img = Image.fromarray(image)
    img.save(outputfile)

    if return_image:
         return image
