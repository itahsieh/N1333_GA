# beam size(")      source1(4A1)/Jy     source2 (4A2)/Jy
# number of wavelength
freq = [22.5e9,43.3e9,100e9,230e9,355e9]
c = 299792458. # m/s
wavelen = [c/x for x in freq]
luminance = []
# 22.5 GHz 
luminance.append( [ [ 0.4,     0.33,    1.554e-3,    0.225e-3] ] )

# 43.3 GHz
luminance.append( [ [0.58,    0.52,    5.57e-3,     0.39e-3] ] ) 

# 100 GHz
luminance.append( [ [1.,       1.,       144e-3,      40e-3,] ] )

# 230 GHz	
luminance.append( [ [1.,       1.,       1230e-3,     369e-3,] ] )

# 355 GHz   
luminance.append(
    [   [0.25,    0.25,    0.388,       0.286],
        [0.4,     0.4,     0.686,       0.393],
        [0.6,     0.6,     0.988,       0.486],
        [0.8,     0.8,     1.194,       0.541],
        [1.0,     1.0,     1.333,       0.579],
        [1.25,    1.25,    1.455,       0.621] ]
)
    

# distance
distance = '235pc'
cell='0.3asec'
npix=32
subres="[['-10asec','-10asec','10asec','10asec',2]]"
