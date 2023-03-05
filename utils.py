import numpy as np

def dms2dd(degrees, minutes, seconds, direction):
    """
    Convert degrees, minutes, seconds to decimal degrees
    """
    dd = float(degrees) + float(minutes)/60 + float(seconds)/(60*60);
    if direction == 'S' or direction == 'W':
        dd *= -1
    return dd;

def dd2dms(deg):
    """
    Convert decimal degrees to degrees, minutes,seconds
    """
    d = np.trunc(deg)
    md = abs(deg - d) * 60
    m = np.trunc(md)
    sd = (md - m) * 60
    return [d, m, sd]


def Marc(bf0, n, PHI0, PHI):
    """
    Compute meridional arc
    
    Input:
        - bf0  : ellipsoid semi major axis multiplied by central meridian 
                 scale factor [meters]; 
        - n    : computed from a, b and f0; 
        - PHI0 : lat of false origin [radians]
        - PHI  : initial or final latitude of point [radians]
    """

    Marc = bf0 * (((1 + n + ((5 / 4) * (n**2)) + ((5 / 4) * (n**3))) * (PHI - PHI0))-
    (((3 * n) + (3 * (n**2)) + ((21 / 8) * (n**3))) * (np.sin(PHI - PHI0)) * (np.cos(PHI + PHI0)))+ 
    ((((15 / 8) * (n**2)) + ((15 / 8) * (n**3))) * (np.sin(2 * (PHI - PHI0))) * (np.cos(2 * (PHI + PHI0))))-
    (((35 / 24) * (n**3)) * (np.sin(3 * (PHI - PHI0))) * (np.cos(3 * (PHI + PHI0)))))

    return Marc

def InitialLat(North, n0, afo, PHI0, n, bfo):
    """
    Compute initial value for Latitude (PHI) [radians]

    Input: 
        - North : northing of point [meters]
        - n0    : northing of false origin (n0) [meters];
        - af0   : ellipsoid semi minor axis multiplied by central meridian 
                 scale factor [meters]; 
        - PHI0  : latitude of false origin (PHI0) [radians]
        - n     : computed from a, b and f0 
        - bf0   : ellipsoid semi major axis multiplied by central meridian 
                 scale factor [meters]; 
    """
    # First PHI value (PHI1)
    PHI1 = ((North - n0) / afo) + PHI0
    
    # Calculate M
    M = Marc(bfo, n, PHI0, PHI1)
    
    # Calculate new PHI value (PHI2)
    PHI2 = ((North - n0 - M) / afo) + PHI1
    
    # Iterate to get final value for InitialLat
    while np.abs(North - n0 - M) > 0.00001:
        PHI2 = ((North - n0 - M) / afo) + PHI1
        M = Marc(bfo, n, PHI0, PHI2)
        PHI1 = PHI2
    # end loop

    InitialLat = PHI2
        
    return InitialLat

def E_N_to_t_minus_T(AtEast, AtNorth, ToEast, ToNorth, a, b, e0, n0, f0, PHI0):
    """
    Compute (t-T) correction in decimal degrees at point (AtEast, AtNorth) 
    to point (ToEast,ToNorth)

    Input:
        - AtEast   : Eastings of point where (t-T) is being
                     being computed [meters], 
        - AtNorth  : Northings of point where (t-T) is 
                     being computed [meters], 
        - ToEast   : Eastings of point at other end of line 
                     to which (t-T) is being computed [meters]
        - ToNorth  : Northings of point at other end of line 
                     to which (t-T) is being computed [meters]
        - a,b      : ellipsoid axis dimensions
        - e0,n0    : easting & northing of true origin in [meters]
        - f0       : central meridian scale factor (f0)
        - PHI0     : latitude of central meridian (PHI0) [decimal degrees]
    """
    # Convert angle measures to radians
    RadPHI0 = PHI0 * (np.pi / 180)
    
    # Compute af0, bf0, e squared (e2), n and Nm (Northing of mid point)
    af0 = a * f0
    bf0 = b * f0
    e2 = ((af0**2) - (bf0**2)) / (af0**2)
    n = (af0 - bf0) / (af0 + bf0)
    Nm = (AtNorth + ToNorth) / 2
    
    # Compute initial value for latitude (PHI) in radians
    PHId = InitialLat(Nm, n0, af0, RadPHI0, n, bf0)
    
    # Compute nu, rho and eta2 using value for PHId
    nu = af0 / (np.sqrt(1 - (e2 * ((np.sin(PHId))**2))))
    rho = (nu * (1 - e2)) / (1 - (e2 * (np.sin(PHId))**2))
    
    # Compute (t-T)
    XXIII = 1 / (6 * nu * rho)
    
    E_N_to_t_minus_T = (180 / np.pi) * ((2 * (AtEast - e0)) + (ToEast - e0)) * (AtNorth - ToNorth) * XXIII
    
    
    return E_N_to_t_minus_T


def E_N_to_C(East, North, a, b, e0, n0, f0, PHI0):
    """
    Compute convergence (in decimal degrees) from easting and northing
    
    Input:
        - East: Eastings [meters]
        - North: Northingsin [meters]
        - a : ellipsoid axis dimensions (a & b) [meters]
        - e0, n0, easting (e0) and northing (n0) of true origin in [meters]
        - f0 : central meridian scale factor
        - PHI0 : latitude of central meridian [decimal degrees]

    """

    # Convert angle measures to radians
    RadPHI0 = PHI0 * (np.pi / 180)
        
    # Compute af0, bf0, e squared (e2), n and Et
    af0 = a * f0
    bf0 = b * f0
    e2 = ((af0**2) - (bf0**2)) / (af0**2)
    n = (af0 - bf0) / (af0 + bf0)
    Et = East - e0
    
    # Compute initial value for latitude (PHI) in radians
    PHId = InitialLat(North, n0, af0, RadPHI0, n, bf0)
    
    # Compute nu, rho and eta2 using value for PHId
    nu = af0 / (np.sqrt(1 - (e2 * ((np.sin(PHId))**2))))
    rho = (nu * (1 - e2)) / (1 - (e2 * (np.sin(PHId))**2))
    eta2 = (nu / rho) - 1

    # Compute Convergence
    XVI = (np.tan(PHId)) / nu
    XVII = ((np.tan(PHId)) / (3 * (nu**3))) * (1 + ((np.tan(PHId))**2) - eta2 - (2 * (eta2**2)))
    XVIII = ((np.tan(PHId)) / (15 * (nu**5))) * (2 + (5 * ((np.tan(PHId))**2)) + (3 * ((np.tan(PHId))**4)))
    
    E_N_to_C = (180 / np.pi) * ((Et * XVI) - ((Et**3) * XVII) + ((Et**5) * XVIII))

    return E_N_to_C


def TrueAzimuth(AtEast, AtNorth, ToEast, ToNorth, a, b, e0, n0, f0, PHI0):
    """
    Compute true azimuth in decimal degrees at point (AtEast, AtNorth) to point (ToEast,ToNorth)
    Input:
       - AtEast   : Eastings of point where (t-T) is being
                     being computed [meters], 
        - AtNorth  : Northings of point where (t-T) is 
                     being computed [meters], 
        - ToEast   : Eastings of point at other end of line 
                     to which (t-T) is being computed [meters]
        - ToNorth  : Northings of point at other end of line 
                     to which (t-T) is being computed [meters]
        - a,b      : ellipsoid axis dimensions
        - e0,n0    : easting & northing of true origin in [meters]
        - f0       : central meridian scale factor (f0)
        - PHI0     : latitude of central meridian (PHI0) [decimal degrees]

    """
    # Compute eastings and northings differences
    Diffe = ToEast - AtEast
    Diffn = ToNorth - AtNorth

    # Compute grid bearing
    if Diffe == 0:
        if Diffn < 0:
            GridBearing = 180
        else:
            GridBearing = 0
        #
    else:
        Ratio = Diffn / Diffe
        GridAngle = (180 / np.pi) * np.arctan(Ratio)
        
        if Diffe > 0:
            GridBearing = 90 - GridAngle
        #
        
        if Diffe < 0:
            GridBearing = 270 - GridAngle
        #
    # End Of ComputeBearing

    # Compute convergence
    Convergence = E_N_to_C(AtEast, AtNorth, a, b, e0, n0, f0, PHI0)
    
    # Compute (t-T) correction
    t_minus_T = E_N_to_t_minus_T(AtEast, AtNorth, ToEast, ToNorth, a, b, e0, n0, f0, PHI0)

    # Compute initial azimuth
    InitAzimuth = GridBearing + Convergence - t_minus_T
    
    # Set TrueAzimuth >=0 and <=360
    if InitAzimuth < 0:
        TrueAzimuth = InitAzimuth + 360
    elif InitAzimuth > 360:
        TrueAzimuth = InitAzimuth - 360
    else:
        TrueAzimuth = InitAzimuth
    #
    return TrueAzimuth

def Lat_Long_to_LSF(PHI, LAM, LAM0, a, b, f0):
    """
    Compute local scale factor from latitude and longitude
    Input:
        PHI: latitude 
        LAM: longitude 
        LAM0: longitude of false origin [decimal degrees]
        a,b: ellipsoid axis dimensions [meters]
        f0: central meridian scale factor
    """
    # Convert angle measures to radians
    RadPHI = PHI * (np.pi / 180)
    RadLAM = LAM * (np.pi  / 180)
    RadLAM0 = LAM0 * (np.pi  / 180)
        
    # Compute af0, bf0 and e squared (e2)
    af0 = a * f0
    bf0 = b * f0
    e2 = ((af0**2) - (bf0**2)) / (af0**2)
    
    # Compute nu, rho, eta2 and p
    nu = af0 / (np.sqrt(1 - (e2 * ((np.sin(RadPHI))**2))))
    rho = (nu * (1 - e2)) / (1 - (e2 * (np.sin(RadPHI))**2))
    eta2 = (nu / rho) - 1
    p = RadLAM - RadLAM0

    # Compute local scale factor
    XIX = ((np.cos(RadPHI)**2) / 2) * (1 + eta2)
    XX = ((np.cos(RadPHI)**4) / 24) * (5 - (4 * ((np.tan(RadPHI))**2)) + (14 * eta2) - (28 * ((np.tan(RadPHI * eta2))**2)))
    
    Lat_Long_to_LSF = f0 * (1 + ((p**2) * XIX) + ((p**4) * XX))

    return Lat_Long_to_LSF

def E_N_to_LSF(East, North, a, b, e0, n0, f0, PHI0):
    """
    Compute local scale factor from from easting and northing
    Input:
        East: Eastings [meters] 
        North: Northings [meters] 
        a,b: ellipsoid axis dimensions (a & b) [meters]
        e0: easting and northing of true origin in [meters]
        f0: central meridian scale factor
        PHI0: latitude of central meridian [decimal degrees]
    """

    # Convert angle measures to radians
    RadPHI0 = PHI0 * (np.pi / 180)
        
    # Compute af0, bf0, e squared (e2), n and Et
    af0 = a * f0
    bf0 = b * f0
    e2 = ((af0**2) - (bf0**2)) / (af0**2)
    n = (af0 - bf0) / (af0 + bf0)
    Et = East - e0
    
    # Compute initial value for latitude (PHI) in radians
    PHId = InitialLat(North, n0, af0, RadPHI0, n, bf0)
    
    # Compute nu, rho and eta2 using value for PHId
    nu = af0 / (np.sqrt(1 - (e2 * ((np.sin(PHId))**2))))
    rho = (nu * (1 - e2)) / (1 - (e2 * (np.sin(PHId))**2))
    eta2 = (nu / rho) - 1

    # Compute local scale factor
    XXI = 1 / (2 * rho * nu)
    XXII = (1 + (4 * eta2)) / (24 * (rho**2) * (nu**2))
    
    E_N_to_LSF = f0 * (1 + ((Et**2) * XXI) + ((Et**4) * XXII))
    
    return E_N_to_LSF
