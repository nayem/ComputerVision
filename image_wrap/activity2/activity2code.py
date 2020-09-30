from numpy import loadtxt
import numpy as np
import scipy.misc as smp
from PIL import Image

def part2(data):
    data = np.asarray(data)
    
    data1 = np.ones((data.shape[0], data.shape[1] + 1))
    data1[:,:-1] = data
    
    transformation_matrix = np.matrix('1 0 0 0; 0 1 0 0; 0 0 1 -5')
    roation_matrix = np.matrix('1 0 0 0; 0 0.71 0.71 0; 0 -0.71 0.71 .25')
    d = 2.5 # 1/2.5 = 0.4
    projection_matrix = np.matrix('1 0 0 0; 0 1 0 0; 0 0 -0.4 0')
    
    new_data = (transformation_matrix * data1.transpose()).transpose()
    xmax = new_data[:, [0]].max()
    ymax = new_data[:, [1]].max()
    print xmax, ymax    
    
    data1[:,:-1] = new_data
    new_data = (projection_matrix * data1.transpose()).transpose()
    xmax = new_data[:, [0]].max()
    ymax = new_data[:, [1]].max()
    print xmax, ymax    
    
    right_eye = np.divide(new_data[:, [0,1]] , new_data[:, [2]]).astype(int)
    xmax = right_eye[:, [0]].max()
    ymax = right_eye[:, [1]].max()
    print xmax, ymax
    
    image1 = np.zeros((ymax+5, xmax+5))
    
    image1[right_eye[:,1], right_eye[:,0]] = 255
    
    img = Image.fromarray(image1, 'L')
    img.show()    

def part3(data):
    data = np.asarray(data)
    
    data1 = np.ones((data.shape[0], data.shape[1] + 1))
    data1[:,:-1] = data
    
    #transformation_matrix = np.matrix('1 0 0 0; 0 1 0 0; 0 0 1 -5')
    transformation_matrix = np.matrix('0.87 0.5 0 0; -0.5 0.87 0 0; 0 0 1 -5') #30 degree rotation
    d = 2.5 # 1/2.5 = 0.4
    projection_matrix = np.matrix('1 0 0 0; 0 1 0 0; 0 0 -0.4 0')
    
    new_data = (transformation_matrix * data1.transpose()).transpose()
    xmax = new_data[:, [0]].max()
    ymax = new_data[:, [1]].max()
    print xmax, ymax    
    
    data1[:,:-1] = new_data
    new_data = (projection_matrix * data1.transpose()).transpose()
    xmax = new_data[:, [0]].max()
    ymax = new_data[:, [1]].max()
    print xmax, ymax    
    
    right_eye = np.divide(new_data[:, [0,1]] , new_data[:, [2]]).astype(int)
    xmax = right_eye[:, [0]].max()
    ymax = right_eye[:, [1]].max()
    print xmax, ymax
    
    image1 = np.zeros((ymax+5, xmax+5))
    
    image1[right_eye[:,1], right_eye[:,0]] = 255
    
    img = Image.fromarray(image1, 'L')
    img.show() 
    
def getEyePlane(data, transformation_matrix, d):    
    
    data = np.asarray(data)    
    data1 = np.ones((data.shape[0], data.shape[1] + 1))
    data1[:,:-1] = data
    #d = 2.5 # 1/2.5 = 0.4
    #projection_matrix = np.matrix('1 0 0 0; 0 1 0 0; 0 0 -0.4 0')
    #print projection_matrix
    projection_matrix = np.matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, -(1.0/d), 0]])
    print projection_matrix
    new_data = (transformation_matrix * data1.transpose()).transpose()
    xmax = new_data[:, [0]].max()
    ymax = new_data[:, [1]].max()
    print xmax, ymax    
    
    data1[:,:-1] = new_data
    new_data = (projection_matrix * data1.transpose()).transpose()
    xmax = new_data[:, [0]].max()
    ymax = new_data[:, [1]].max()
    print xmax, ymax    
    
    right_eye = np.divide(new_data[:, [0,1]] , new_data[:, [2]]).astype(int)
    
    return right_eye

def part4(data):   
    
    transformation_matrix = np.matrix('0 1 0 0; -1 0 0 0; 0 0 1 -50')
    right_eye = abs(getEyePlane(data, transformation_matrix, 40))
    xmax0 = right_eye[:, [0]].max()
    ymax0 = right_eye[:, [1]].max()
    print xmax0, ymax0
    
    transformation_matrix = np.matrix('0 1 0 30; -1 0 0 0; 0 0 1 -50')
    left_eye = abs(getEyePlane(data, transformation_matrix, 40))
    xmax1 = left_eye[:, [0]].max()
    ymax1 = left_eye[:, [1]].max()
    print xmax1, ymax1   
    
    xmax = max(xmax0, xmax1)
    ymax = max(ymax0,ymax1)
    
    image0 = np.zeros((ymax+5, xmax+5))    
    image0[right_eye[:,1], right_eye[:,0]] = 255
    
    image1 = np.zeros((ymax+5, xmax+5))    
    image1[left_eye[:,1], left_eye[:,0]] = 255    
    
    image = np.zeros((ymax+5, xmax+5, 3), dtype=np.uint8)
    image[:, :, 0] = image0;
    image[:, :, 1] = image1;
    image[:, :, 2] = image1;
    img = Image.fromarray(image, 'RGB')
    #img.show()
    img.save('my.png') 
        
# part 1
ins = open( "b657-wars.txt", "r" )
data = []
for line in ins:
    number_strings = line.split() # Split the line on runs of whitespace
    numbers = [int(n) for n in number_strings] # Convert to integers
    data.append(numbers) # Add the "row" to your list.
ins.close()
#print data

#part 2
#part2(data)

#part 3
#part3(data)

#part 4
part4(data)