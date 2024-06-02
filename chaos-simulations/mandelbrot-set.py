from PIL import Image

#This code is based on code from Math Adventures with python. p.132

def mandelbrot(z,num):
    zIter = z
    for count in range(num+1):
        if abs(zIter) > 2:
            return count
        zIter = zIter**2 + z
    return num


halfWidth = 300
xPos = 400
yPos = 300
zoom = halfWidth/2


img = Image.new( 'HSV', (halfWidth*2,halfWidth*2), "white") 
pixels = img.load() 
for i in range(img.size[0]):    
    for j in range(img.size[1]):
        zij = complex(i-xPos,j-yPos)/(zoom)
        m =mandelbrot(zij,100)
        pixels[i,j] = (225-5*m, 225, 225) 

img.show()
