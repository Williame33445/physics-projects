from abc import abstractmethod, ABC
from particle import *

"""if the tree increases to the max depth then the objects tree becomes a non breaching node,
which is basically a list of the particles inside it"""
def treeNode(midPoint,halfWidth,maxDepth=15,theta=0.01):
    if maxDepth > 1:
        return BranchingNode(midPoint,halfWidth,maxDepth,theta)
    else:
        return NonBranchingNode(midPoint,halfWidth,theta)

"""the parents class that has abstract methods that need to be defined in the subclasses"""
class AbstractNode(ABC):
    def __init__(self,midPoint,halfWidth,theta):
        self.midPoint = midPoint
        self.halfWidth = halfWidth
        self.theta = theta
        self.combinedParticle = None

    @abstractmethod
    def addParticle(self,newParticle):
        pass

    @abstractmethod
    def findMassDistribution(self):
        pass
    """finds the acceleration of the particle through recursion of calculate net acceleration until its 
    below the theta value or a leaf node"""
    def netAccelerationOf(self,targetParticle):
        r = self.combinedParticle.pos.distanceTo(targetParticle.pos)
        d = self.halfWidth * 2
        if (d < self.theta * r):
            return targetParticle.accelerationTowards(self.combinedParticle)
        else:
            return self.calculateNetAcceleration(targetParticle)

    @abstractmethod
    def calculateNetAcceleration(self,targetParticle):
        pass


# A node that can create children when particles are added to it
class BranchingNode(AbstractNode):
    def __init__(self,midPoint,halfWidth,maxDepth,theta):
        AbstractNode.__init__(self,midPoint,halfWidth,theta)
        self.maxDepth = maxDepth
        self.childNodes = [None,None,None,None]
        self.particleCount = 0
    def children(self):
        return filter(None,self.childNodes)
    def addParticle(self,newParticle):
        #if there are no particles in this node put the particle in here
        if self.particleCount == 0:
            self.combinedParticle = newParticle
        else:
            """if there are other particles in this node then find the segment
            that the new particle goes in and add to node"""
            #add particle to child node
            self.addToCorrectChild(newParticle)
            """this recusively runs until we get to a leaf node(bottom node) and extends it
            until one bellow its own node"""
            if self.combinedParticle != None:
                #does the same proccess for the existing particle until both particles have their own node
                self.addToCorrectChild(self.combinedParticle)
                #wipe the particle as we now have multiple particles in the same node
                self.combinedParticle = None

        self.particleCount += 1

    def addToCorrectChild(self,particle):
        #find quadrent that the particle should go in
        childIndex = self.quadrantNumberFor(particle.pos)
        """if the child isn't occupied it returns none and calculates the childs midpoint and halfWidth
        otherwise it returns the values of the particles that occupy it"""
        child = self.getOrCreateChildAt(childIndex)
        """recursion then happens again and the entire proccess continues to happen until the particles
        each have there own node"""
        child.addParticle(particle)

    # Function that finds out what segment position is in: 
    #     0 - top left, 1 - top right, 2 bottom left, 3 - bottom right
    def quadrantNumberFor(self,position):
        midPoint = self.midPoint
        if (midPoint.x>position.x)and(midPoint.y>position.y):
            return 0
        elif (midPoint.x<=position.x)and(midPoint.y>position.y):
            return 1
        elif (midPoint.x>position.x)and(midPoint.y<=position.y):
            return 2
        else:
            return 3

    def getOrCreateChildAt(self, childIndex):
        #we define (0,0) at the centre of the screen
        #if the child index is not occupied
        if self.childNodes[childIndex] == None:
            #calculate the width of the child
            childHalfWidth = self.halfWidth / 2.0
            """calculate the new child's midpoint by adding or subtracting the
            appropriate offsets"""
            #are defined due to the way segments were defined
            deltaX = [-childHalfWidth, +childHalfWidth, -childHalfWidth, +childHalfWidth]
            deltaY = [-childHalfWidth, -childHalfWidth, +childHalfWidth, +childHalfWidth]
            #calculate the offset
            offset = Vector(deltaX[childIndex],deltaY[childIndex])
            childMidpoint = self.midPoint.plus(offset)
            #define the child
            self.childNodes[childIndex] = treeNode(childMidpoint,childHalfWidth,self.maxDepth-1,self.theta)
        return self.childNodes[childIndex]

    def findMassDistribution(self):
        def childCentresOfMass():
            result = []
            for c in self.children():
                result.append(c.combinedParticle)
            return result  
        #if the object is a leaf then no calculations are required
        if self.combinedParticle != None:
            return
        
        for c in self.children():
            #recursively runs find mass distribtion through the children of the nodes to find mass and center of mass
            c.findMassDistribution()
  
        self.combinedParticle = combinedParticle(childCentresOfMass())

    def calculateNetAcceleration(self,targetParticle):
        # If this is a leaf node (no children)..
        if self.particleCount == 1:
            return targetParticle.accelerationTowards(self.combinedParticle)
        else:
            netAcceleration = zeroVector()
            for c in self.children():
                # Make a mutually recursive call for each child
                netAcceleration.translate(c.netAccelerationOf(targetParticle))
            return netAcceleration


# A node will not create children when particles are added to it but will hold them in a list
class NonBranchingNode(AbstractNode):
    def __init__(self,midPoint,halfWidth,theta):
        AbstractNode.__init__(self,midPoint,halfWidth,theta)
        self.particles = []

    def addParticle(self,newParticle):
        self.particles.append(newParticle)

    def findMassDistribution(self):
        #finds the center of mass averages the parameters for the nodes
        self.combinedParticle = combinedParticle(self.particles)

    #if this gets called then its bellow the theta value and the particles are calculated exactly
    def calculateNetAcceleration(self,targetParticle):
        netAcceleration = zeroVector()
        for p in self.particles:
            netAcceleration.translate(targetParticle.accelerationTowards(p))
        return netAcceleration
