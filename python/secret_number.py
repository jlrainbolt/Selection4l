def generate_secret_number():
    import struct
    import numpy as np
    r = np.random.random()
    a = 0.0
    if r > 0.5 :
        a = np.random.normal(0.05, 0.025)
    else : 
        a = np.random.normal(-0.05, 0.025)
    newfile = open("secret_number.txt", "wb")
    newfile.write(struct.pack('f', a))
    newfile.close()

def read_secret_number():
    newfile = open("../data/secret_number.txt", "rb")
    a = newfile.read()
    import struct
    a = struct.unpack('f', a)
    return a[0]
