from particle import zeroVector
from tree import treeNode

#Class that deals with running the simulation
class Simulator:

    def __init__(self,particles,halfWidth,theta,maxDepth):
        #list of all particles 
        self.particles = particles
        #half width of the range that we are thinking about
        self.halfWidth = halfWidth
        self.theta = theta
        self.maxDepth = maxDepth

    #creates the tree
    def buildTree(self):
        self.rootNode = treeNode(zeroVector(),self.halfWidth,self.maxDepth,self.theta)
        self.addParticles()
        self.rootNode.findMassDistribution()

    #adds all particles
    def addParticles(self):
        for p in self.particles:
            self.rootNode.addParticle(p)
    #finds acceleration of a single particle
    def netAccelerationOf(self,particle):
        return self.rootNode.netAccelerationOf(particle)

    #finds the accelerations of all the particles and then applies them
    def applyAcceleration(self,timeElapsed):
        def allAccelerations():
            result = []
            for p in self.particles:
                result.append(self.netAccelerationOf(p))
            return result

        accelerations = allAccelerations()
        for x in range(len(self.particles)):
            acceleration = accelerations[x]
            deltaVelocity = acceleration.times(timeElapsed)
            self.particles[x].velocity.translate(deltaVelocity)
        
    #applies the velocity to change the position
    def applyVelocity(self,timeElapsed):
        for p in self.particles:
            deltaPosition = p.velocity.times(timeElapsed)
            p.pos.translate(deltaPosition)

    #cycles through rebuilding the tree and applying the velocities
    def tick(self,timeElapsed):
        self.buildTree()
        self.applyVelocity(timeElapsed)
        self.applyAcceleration(timeElapsed)

#Class the controls how long the program runs for
class SimulationParams:
    def __init__(self,halfWidth,tickPeriod,totalDuration,theta,maxDepth):
        self.halfWidth = halfWidth
        self.tickPeriod = tickPeriod
        self.totalDuration = totalDuration
        self.theta = theta
        self.maxDepth = maxDepth

    def numberOfCycles(self):
        return int(self.totalDuration/self.tickPeriod)

    def width(self):
        return 2.0 * self.halfWidth
