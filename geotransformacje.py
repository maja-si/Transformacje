from math import sin, cos, sqrt, tan, atan, atan2, degrees, radians, pi
from numpy import array
class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2


    
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, defoulf 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
            

    def geod2XYZ(self, fi, lam, h):
        
        N = self.a/sqrt(1-self.ecc2*sin(fi)**2)
        X = (N + h) * cos(fi) * cos(lam)
        Y = (N + h) * cos(fi) * sin(lam)
        Z = (N*(1-self.ecc2) + h) * sin(fi)
        return(X, Y, Z)

    def u1992(self, fi, lam):
        m_0 = 0.9993
        N = self.a/sqrt(1-self.ecc2*sin(fi)**2)
        t = tan(fi)
        e_2 = self.ecc2/(1-self.ecc2)
        n2 = e_2 * (cos(fi))**2
    
        lam_0 = radians(19)
        l = lam - lam_0
    
        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)   
        A2 = (3/8) * (self.ecc2 + ((self.ecc2**2)/4) + ((15 * (self.ecc2**3))/128))
        A4 = (15/256) * (self.ecc2**2 + ((3*(self.ecc2**3))/4))
        A6 = (35 * (self.ecc2**3))/3072 
    

        sig = self.a * ((A0*fi) - (A2*sin(2*fi)) + (A4*sin(4*fi)) - (A6*sin(6*fi)))
    
        x = sig + ((l**2)/2) * N *sin(fi) * cos(fi) * (1 + ((l**2)/12) * ((cos(fi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((l**4)/360) * ((cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*cos(fi)) * (1 + ((((l**2)/6) * (cos(fi))**2) * (1-t**2+n2)) +  (((l**4)/(120)) * (cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x92 = round(m_0*x - 5300000, 3)
        y92 = round(m_0*y + 500000, 3)
        return x92, y92 

    def u2000(self, fi, lam):
        m_0 = 0.999923
        N = self.a/sqrt(1-self.ecc2*sin(fi)**2)
        t = tan(fi)
        e_2 = self.ecc2/(1-self.ecc2)
        n2 = e_2 * (cos(fi))**2
    
        lam = degrees(lam)
        if lam>13.5 and lam <16.5:
            s = 5
            lam0 = 15
        elif lam>16.5 and lam <19.5:
            s = 6
            lam0 = 18
        elif lam>19.5 and lam <22.5:
            s = 7
            lam0 = 21
        elif lam>22.5 and lam <25.5:
            s = 8
            lam0 = 24
        
        lam = radians(lam)
        lam0 = radians(lam0)
        l = lam - lam0

        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)   
        A2 = (3/8) * (self.ecc2 + ((self.ecc2**2)/4) + ((15 * (self.ecc2**3))/128))
        A4 = (15/256) * (self.ecc2**2 + ((3*(self.ecc2**3))/4))
        A6 = (35 * (self.ecc2**3))/3072 
    

        sig = self.a * ((A0*fi) - (A2*sin(2*fi)) + (A4*sin(4*fi)) - (A6*sin(6*fi)))
    
        x = sig + ((l**2)/2) * N *sin(fi) * cos(fi) * (1 + ((l**2)/12) * ((cos(fi))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((l**4)/360) * ((cos(fi))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        y = l * (N*cos(fi)) * (1 + ((((l**2)/6) * (cos(fi))**2) * (1-t**2+n2)) +  (((l**4)/(120)) * (cos(fi)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
        
        x00 = round(m_0 * x, 3) 
        y00 = round(m_0 * y + (s*1000000) + 500000, 3)
     
        return(x00, y00)


    def distance2D (self, xA, yA, xB, yB):
        dX = xB - xA
        dY = yB - yA
        dist = sqrt(dX**2+dY**2)
        return dist
    
    def distance3D(self, xA, yA, xB, yB, zA, zB):
        dX = xB - xA
        dY = yB - yA
        dZ = zB - zA
        dist = sqrt(dX**2+dY**2+dZ**2)
        return dist

    
    def topocentryczne(self, fi, lam, h):
        N = self.a/sqrt(1-self.ecc2*sin(fi)**2)
        X = (N + h) * cos(fi) * cos(lam)
        Y = (N + h) * cos(fi) * sin(lam)
        Z = (N*(1-self.ecc2) + h) * sin(fi)
        return(X, Y, Z)
        R = array([[-sin(fi) * cos(lam), -sin(lam), cos(fi)*cos(lam)],
                   [-sin(fi)*sin(lam), cos(lam), cos(fi)*sin(lam)],
                   [cos(fi), 0, sin(fi)]])
        neu = (R.transpose()) @ array([X, Y, Z])
        return neu
    def azimuth(self, xA, yA, xB, yB):
        dX = xB - xA
        dY = yB - yA
        Az_rad = atan2(dY, dX)
        Az_deg = degrees(Az_rad)
        if Az_deg < 0:
            Az_deg = Az_deg+ 360
        return Az_deg
        
if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "wgs84")
    # dane XYZ geocentryczne
    X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    phi, lam, h = geo.xyz2plh(X, Y, Z)
    print("fi=", phi, "lam=", lam, "h=", h)

