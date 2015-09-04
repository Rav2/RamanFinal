from pylab import *
import matplotlib.colors


def middle_point(data):
        av_norm=0.5
        av_data=average(data)
        length=0
        for i in data:
            for j in i:
                length+=1
        if length !=0:
            av_norm=av_data/length
            return av_norm
        else:
            return None

def minimum(data):
    minimal=data[0][0]
    for i in data:
        for j in i:
            if j<minimal:
                minimal=j
    return minimal

def maximum(data):
    maximal=data[0][0]
    for i in data:
        for j in i:
            if j>maximal:
                maximal=j
    return maximal

def make_map(data,vmin=0.0,vmax=1.0,av=0.5,):
    data_array=np.array(data)
    #USTAWIENIE KOLOROWANIA
    cdict = {'red': ((0.0, 0.0, 0.0),
                    (av, 1.0, 1.0),
                     (1.0, 0.0, 0.0)),
             'green': ((0.0, 0.0, av),
                       (av, 1.0, 1.0),
                       (1.0, av, 0.0)),
             'blue': ((0.0, 0.0, 0.0),
                       (av, 1.0, 1.0),
                       (1.0, 0.0, 0.0))}
    #UTWORZENIE MAPY I PLIKU GRAFICZNEGO
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    pcolor(data_array,cmap=my_cmap,vmin=vmin,vmax=vmax)
    colorbar()
    #show()
    #savefig("Map.png")

