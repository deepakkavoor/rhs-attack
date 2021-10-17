import random
# max breadth, height of outer quadrant
# change
b = 2**17
h = 2**17
# b=16
# h=16
xmin = 0
xmax = b
ymin = 0
ymax = h

def get_quadtree(m):
    # return quadtree N with max m nodes; each parent has 4 children
    N = []

    # first node: biggest outer quadrant
    # Ni1 starts from top left and Ni4 ends at bottom left
    # x increases from left->right and y increases from bottom->top
    i=0 
    Ni1 = [xmin, ymax]
    Ni2 = [xmax, ymax]
    Ni3 = [xmax, ymin]
    Ni4 = [xmin, ymin]
    Ni = [Ni1,Ni2,Ni3,Ni4]
    N.append(Ni)

    while len(N)+4 <= m:
        i+=1
        N_prev = N[i-1]
        # choose a random center within N_prev
        x = random.randint(N_prev[0][0], N_prev[1][0])
        y = random.randint(N_prev[3][1], N_prev[0][1])    
        # divide N_prev into 4 quadrants

        # topleft
        Ni1 = N_prev[0]
        Ni2 = [x, N_prev[0][1]]
        Ni3 = [x, y]
        Ni4 = [N_prev[0][0], y]
        Ni = [Ni1,Ni2,Ni3,Ni4]
        N.append(Ni)

        # topright
        Ni1 = [x, N_prev[1][1]]
        Ni2 = N_prev[1]
        Ni3 = [N_prev[1][0], y]
        Ni4 = [x, y]
        Ni = [Ni1,Ni2,Ni3,Ni4]
        N.append(Ni)

        # bottomright
        Ni1 = [x, y]
        Ni2 = [N_prev[2][0], y]
        Ni3 = N_prev[2]
        Ni4 = [x, N_prev[2][1]]
        Ni = [Ni1,Ni2,Ni3,Ni4]
        N.append(Ni)

        # bottomleft
        Ni1 = [N_prev[3][0], y]
        Ni2 = [x,y]
        Ni3 = [x, N_prev[3][1]]
        Ni4 = N_prev[3]
        Ni = [Ni1,Ni2,Ni3,Ni4]
        N.append(Ni)
    return N

if __name__ == '__main__':
    m = random.randint(28, 84)
    N = get_quadtree(m)
    print(len(N))