import tables

def get_stuff(f):

    # get the charge
    charge = f.get_node('/Homogenized_QTot')
    charge = charge.col('value')

    # weight dict and prop matrix
    weight_dict = f.get_node('/JulietWeightDict')
    prop_matrix = [f.get_node(f'/PropagationMatrix{flav}') for flav in ['NuE', 'NuMu', 'NuTau']]

    evts_per_file = weight_dict.col('EventsPerFile')[0]

    return charge, weight_dict, prop_matrix, evts_per_file

def find_contours_2D(x_values,y_values,xbins,weights=None,c1=16,c2=84):   
    """
    Find upper and lower contours and median
    x_values = array, input for hist2d for x axis (typically truth)
    y_values = array, input for hist2d for y axis (typically reconstruction)
    xbins = values for the starting edge of the x bins (output from hist2d)
    c1 = percentage for lower contour bound (16% - 84% means a 68% band, so c1 = 16)
    c2 = percentage for upper contour bound (16% - 84% means a 68% band, so c2=84)
    Returns:
        x = values for xbins, repeated for plotting (i.e. [0,0,1,1,2,2,...]
        y_median = values for y value medians per bin, repeated for plotting (i.e. [40,40,20,20,50,50,...]
        y_lower = values for y value lower limits per bin, repeated for plotting (i.e. [30,30,10,10,20,20,...]
        y_upper = values for y value upper limits per bin, repeated for plotting (i.e. [50,50,40,40,60,60,...]
    """
    if weights is not None:
        import wquantiles as wq
    y_values = numpy.array(y_values)
    indices = numpy.digitize(x_values,xbins)
    r1_save = []
    r2_save = []
    median_save = []
    for i in range(1,len(xbins)):
        mask = indices==i
        if len(y_values[mask])>0:
            if weights is None:
                r1, m, r2 = numpy.percentile(y_values[mask],[c1,50,c2])
            else:
                r1 = wq.quantile(y_values[mask],weights[mask],c1/100.)
                r2 = wq.quantile(y_values[mask],weights[mask],c2/100.)
                m = wq.median(y_values[mask],weights[mask])
        else:
            #print(i,'empty bin')
            r1 = numpy.nan
            m = numpy.nan
            r2 = numpy.nan
        median_save.append(m)
        r1_save.append(r1)
        r2_save.append(r2)
    median = numpy.array(median_save)
    lower = numpy.array(r1_save)
    upper = numpy.array(r2_save)

    # the first return with the [1:] and [:-1] is about locating the bin centers
    return numpy.asarray((xbins[1:] + xbins[:-1])/2), numpy.asarray(median), numpy.asarray(lower), numpy.asarray(upper)