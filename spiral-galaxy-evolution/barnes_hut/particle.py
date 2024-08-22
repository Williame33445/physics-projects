import math 
global G
G = 6.67*10**-11

class Vector:
    def __init__(self,x,y):
        self.x = x
        self.y = y
    def plus(self,delta):
        return Vector(self.x+delta.x, self.y+delta.y)
    def minus(self,delta):
        return Vector(self.x-delta.x, self.y-delta.y)
    def times(self,factor):
        return Vector(self.x*factor, self.y*factor)
    def dividedBy(self,factor):
        return Vector(self.x/factor, self.y/factor)
    def distanceTo(self,pos2):
        return math.sqrt((self.x-pos2.x)**2+(self.y-pos2.y)**2)
    def translate(self,delta):
        self.x += delta.x
        self.y += delta.y
    def __str__(self):
        #allows you to define what str() does for this class
        return "(" + str(self.x) + "," + str(self.y) + ")"

class Particle:
    def __init__(self,mass,pos):
        self.mass = mass
        self.pos = pos
    def __str__(self):
        #allows you to define what str() does for this class
        return "mass=" + str(self.mass) + ",centre=" + str(self.pos)
    def accelerationTowards(self,particle1):
        r = particle1.pos.distanceTo(self.pos)
        if r==0:
            return Vector(0.0, 0.0)
        magnitude = (G)*(particle1.mass)/(r**3)
        #gives vector in the direction of the
        direction = particle1.pos.minus(self.pos)
        return direction.times(magnitude)


#gives the initial state of the system
class KinematicParticle(Particle):
    def __init__(self,mass,pos,velocity):
        Particle.__init__(self,mass,pos)
        self.velocity = velocity
    def __str__(self):
        #allows you to define what str() does for this class
        return Particle.__str__(self) + ",velocity=" + str(self.velocity)

#allows you to create a particle without first having to create a position
def particle(mass,x,y):
    return Particle(mass,Vector(x,y))

def kinematicParticle(mass,x,y,vx,vy):
    return KinematicParticle(mass,Vector(x,y),Vector(vx,vy))


def zeroVector():
    return Vector(0.0, 0.0)
    
def combinedParticle(particles):
    mass = 0
    centre = zeroVector()
    for p in particles:
        mass += p.mass
        centre = centre.plus(p.pos.times(p.mass))
    centre = centre.dividedBy(mass)
    return Particle(mass, centre)

