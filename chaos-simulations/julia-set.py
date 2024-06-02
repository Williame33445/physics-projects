from PIL import Image

#This code is based on code from Math Adventures with python. p.141


def julia(z,num,c):
    zIter = z
    for count in range(num+1):
        if abs(zIter) > 2:
            return count
        zIter = zIter**2 + c
    return num


halfWidth = 300
xPos = 400
yPos = 300
zoom = halfWidth/3


img = Image.new( 'HSV', (halfWidth*2,halfWidth*2), "white") 
pixels = img.load() 
for i in range(img.size[0]):    
    for j in range(img.size[1]):
        zij = complex(i-xPos,j-yPos)/(zoom)
        m =julia(zij,200,complex(-0.4,0.6))
        pixels[i,j] = (225-5*m, 225, 225) 

img.show()
