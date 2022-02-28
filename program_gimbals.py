
from matplotlib import pyplot as plt
import numpy as np
#from mpl_toolkits import mplot3d
from scipy.linalg import expm, norm
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import tkinter as tk


quat1 =[] # It stores the quaternions for inner gimbal (1)
quat2 =[] # It stores the quaternions for middle gimbal (2)
quat3 =[] # It stores the quaternions for outer gimbal (3)

Euler321rad = [] # It stores the Euler angles for the gimbal 1 for 
                 # the convention ZYX or 321 in radians


'''
Each colunm contains the timestep and the quaternion components for the 
simulation.

q0 real 
q1,q2,q3 imaginary

quat1 = [array([Ts_1 q0 q1 q2 q3]     -> first simulation (run)
               [ .   .  .  .  .]
                       .
                       .
                       .
               )                    
         
         array([Ts_2 q0 q1 q2 q3]      -> second simulation (run)
               [ .   .  .  .  .]
                       .
                       .
                       .        
               )                    
               
         ]
'''


def quatMult0123(qa, qb):
    q0 = qa[0]
    q1 = qa[1]
    q2 = qa[2]
    q3 = qa[3]
    q = qb.reshape(4, 1)

    return (
        np.array(
            [
                [q0, -q1, -q2, -q3],
                [q1, q0, -q3, q2],
                [q2, q3, q0, -q1],
                [q3, -q2, q1, q0],
            ]
        ).reshape(4, 4)
        @ q
    )


def skew(w):
    return np.array(
        [[0, -w[2], w[1]], [w[2], 0, -w[0]], [-w[1], w[0], 0]], dtype=object
    )


def inv_quaternion(q):
    return np.array([q[0], -q[1], -q[2], -q[3]]) / (norm(q) ** 2)


def quaternion_kinematics(w, Ts):
    # qp = 0.5 q prod_quat wq
    omega = np.zeros((4, 4))
    omega[0, 1] = -w[0]
    omega[0, 2] = -w[1]
    omega[0, 3] = -w[2]
    omega[1, 0] = w[0]
    omega[1, 2] = w[2]
    omega[1, 3] = -w[1]
    omega[2, 0] = w[1]
    omega[2, 1] = -w[2]
    omega[2, 3] = w[0]
    omega[3, 0] = w[2]
    omega[3, 1] = w[1]
    omega[3, 2] = -w[0]

    return expm(0.5 * omega * Ts)

def simulation(w1, w2, w3, time, Ts, q1, q2, q3, phi0, theta0):
    w1 = w1.reshape(3, 1)
    w2 = w2.reshape(3, 1)
    w3 = w3.reshape(3, 1)

    t = np.arange(0, time, Ts)
    Q1 = np.zeros((4, len(t)))
    Q2 = np.zeros((4, len(t)))
    Q3 = np.zeros((4, len(t)))

    for i in range(len(t)):
        # guimbal 3 -> outer    
        w3I3 = w3 # angular velocity of gimbal 3 in relation to the 
                  #inertial referential, written in system 3
        
        # guimbal 2 -> middle
        theta = theta0 + w2[1] * t[i]
        C32 = np.array(
            [
                [np.cos(theta), 0, -np.sin(theta)],
                [0, 1, 0],
                [np.sin(theta), 0, np.cos(theta)],
            ],
            dtype=object,
        )

        
        w2I2 = w2 + C32 @ w3 # angular velocity of gimbal 2 in relation to the 
                             # inertial referential, written in system 2
        
        
        # guimbal 1 -> inner
        phi = phi0 + w1[0] * t[i]
        C21 = np.array(
            [[1, 0, 0], [0, np.cos(phi), np.sin(phi)], [0, -np.sin(phi), np.cos(phi)]],
            dtype=object,
        )

        w1I1 = w1 + C21 @ w2 + C21 @ C32 @ w3 # angular velocity of gimbal 1 in relation to the 
                                              #inertial referential, written in system 1
        
        # quaternions
        q3 = quaternion_kinematics(w3I3, Ts) @ q3.reshape(4,1)
        q3 = q3 / norm(q3)
        q2 = quaternion_kinematics(w2I2, Ts) @ q2.reshape(4,1)
        q2 = q2 / norm(q2)
        q1 = quaternion_kinematics(w1I1, Ts) @ q1.reshape(4,1)
        q1 = q1 / norm(q1)

        Q1[:, i] = q1.reshape(4,)
        Q2[:, i] = q2.reshape(4,)
        Q3[:, i] = q3.reshape(4,)

    return Q1, Q2, Q3, phi, theta


def rotation(vector, q):
    q = q.reshape(4, 1)
    v = np.block([0, vector]).reshape(4, 1)
    qinv = inv_quaternion(q).reshape(4, 1)
    r = quatMult0123(q, v)
    r = quatMult0123(r, qinv)

    return r[1:]

def quaternion2Euler(q):
    # 321 sequence Inertial to Body frame
    # q0,q1,q2,q3 convention
    
    euler = np.zeros((q.shape[0],3));
    q0 = q[:,0];
    q1 = q[:,1];
    q2 = q[:,2];
    q3 = q[:,3];
    euler[:,2] = np.arctan2( 2 * q1 * q2 + 2 * q0 * q3 , 2 * q0 * q0 + 2 * q1 * q1 - 1);
    euler[:,1] = np.arcsin(-2*q1*q3 + 2*q0*q2);
    euler[:,0] = np.arctan2( 2 * q2 * q3 + 2 * q0 * q1 , 2 * q0 * q0 + 2 * q3 * q3 - 1);
    return euler

class Arrow3D(FancyArrowPatch):

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs
    
    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)



class gimbals:

	def __init__(self):
		
	    # each row of the self.frames represents a point
	    self.framex = np.array([[1, 0, -1], [1, 0, 1], [-1, 0, 1], [-1, 0, -1], [1, 0, -1]])
	    self.framez = 2 * np.array([[0, -1, -1], [0, -1, 1], [0, 1, 1], [0, 1, -1], [0, -1, -1]])
	    self.framey = 1.5 * np.array([[1, 1, 0], [-1, 1, 0], [-1, -1, 0], [1, -1, 0], [1, 1, 0]])

	    self.connect_zy = np.array(
	        [
	            [(self.framez[1, 0] + self.framez[0, 0]) / 2, (self.framez[1, 1] + self.framez[0, 1]) / 2, 0],
	            [(self.framey[2, 0] + self.framey[3, 0]) / 2, (self.framey[2, 1] + self.framey[3, 1]) / 2, 0],
	        ]
	    )
	    self.connect_zy1 = np.array(
	        [
	            [(self.framez[2, 0] + self.framez[3, 0]) / 2, (self.framez[2, 1] + self.framez[3, 1]) / 2, 0],
	            [(self.framey[1, 0] + self.framey[0, 0]) / 2, (self.framey[1, 1] + self.framey[0, 1]) / 2, 0],
	        ]
	    )

	    self.connect_zI = np.array(
	        [
	            [0, 0, (self.framez[1, 2] + self.framez[2, 2]) / 2],
	            [0, 0, self.framez[1, 2] + 0.5],
	        ]
	    )

	    self.connect_zI1 = np.array(
	        [
	            [0, 0, (self.framez[0, 2] + self.framez[3, 2]) / 2],
	            [0, 0, self.framez[0, 2] - 0.5],
	        ]
	    )

	    self.connect_xy = np.array(
	        [
	            [
	                (self.framex[0, 0] + self.framex[1, 0]) / 2,
	                (self.framex[0, 1] + self.framex[1, 1]) / 2,
	                (self.framex[0, 2] + self.framex[1, 2]) / 2,
	            ],
	            [
	                (self.framey[3, 0] + self.framey[0, 0]) / 2,
	                (self.framey[3, 1] + self.framey[0, 1]) / 2,
	                (self.framey[3, 2] + self.framey[0, 2]) / 2,
	            ],
	        ]
	    )

	    self.connect_xy1 = np.array(
	        [
	            [
	                (self.framex[2, 0] + self.framex[3, 0]) / 2,
	                (self.framex[2, 1] + self.framex[3, 1]) / 2,
	                (self.framex[2, 2] + self.framex[3, 2]) / 2,
	            ],
	            [
	                (self.framey[1, 0] + self.framey[2, 0]) / 2,
	                (self.framey[1, 1] + self.framey[2, 1]) / 2,
	                (self.framey[1, 2] + self.framey[2, 2]) / 2,
	            ],
	        ]
	    )

	    self.phi0 = 0
	    self.theta0 = 0
	    a1 = 0
	    a2 = 0
	    a3 = 0
	    n1 = np.array([1, 0, 0])
	    n2 = np.array([0, 1, 0])
	    n3 = np.array([0, 0, 1])

	    self.q3_ = np.block([np.cos(a3 / 2), n3 * np.sin(a3 / 2)])
	    self.q2_ = np.block([np.cos(a2 / 2), n2 * np.sin(a2 / 2)])
	    self.q1_ = np.block([np.cos(a1 / 2), n1 * np.sin(a1 / 2)])
        
        
	    self.body_frame = 0.5*np.array([[1,0,0],[0,0,0],[0,1,0],[0,0,0],[0,0,1]]) 
        
	   
	    self.fig = plt.figure()
	    self.ax = self.fig.add_subplot(111, projection="3d")
       
	   
	    


	def animation(self, phi, theta, psi, time, Ts):
	    
	    w1 = np.array([phi, 0, 0]) / time
	    w2 = np.array([0, theta, 0]) / time
	    w3 = np.array([0, 0, psi]) / time

	    q1, q2, q3, self.phi0, self.theta0 = simulation(w1, w2, w3, time, Ts, self.q1_, self.q2_, self.q3_, self.phi0, self.theta0)

	    framex_ = np.zeros((self.framex.shape[0],self.framex.shape[1],q1.shape[1]))
	    framey_ = np.zeros((self.framex.shape[0],self.framex.shape[1],q1.shape[1]))
	    framez_ = np.zeros((self.framex.shape[0],self.framex.shape[1],q1.shape[1]))
	    connect_zy_ = np.zeros((self.connect_zy.shape[0],self.connect_zy.shape[1],q1.shape[1]))
	    connect_zy1_ = np.zeros((self.connect_zy.shape[0],self.connect_zy.shape[1],q1.shape[1]))
	    connect_xy_ = np.zeros((self.connect_xy.shape[0],self.connect_xy.shape[1],q1.shape[1]))
	    connect_xy1_ = np.zeros((self.connect_xy.shape[0],self.connect_xy.shape[1],q1.shape[1]))
	    body_ = np.zeros((self.body_frame.shape[0],self.body_frame.shape[1],q1.shape[1]))

	   
	    arrow = np.array([0.7,0,0])

	    for i in range(q1.shape[1]):
	        
	        for j in range(len(self.framex)):
	            framex_[j, :, i] = rotation(self.framex[j, :], q1[:, i]).reshape(
	                3,
	            )

	            framey_[j, :, i] = rotation(self.framey[j, :], q2[:, i]).reshape(
	                3,
	            )
	            framez_[j, :, i] = rotation(self.framez[j, :], q3[:, i]).reshape(
	                3,
	            )
	        for k in range(len(self.connect_xy)):
	            connect_zy_[k, :, i] = rotation(self.connect_zy[k, :], q3[:, i]).reshape(
	                3,
	            )
	            connect_zy1_[k, :, i] = rotation(self.connect_zy1[k, :], q3[:, i]).reshape(
	                3,
	            )
	            connect_xy_[k, :, i] = rotation(self.connect_xy[k, :], q2[:, i]).reshape(
	                3,
	            )
	            connect_xy1_[k, :, i] = rotation(self.connect_xy1[k, :], q2[:, i]).reshape(
	                3,
	            )
	        for j in range(len(self.body_frame)):
	            body_[j, :, i] = rotation(self.body_frame[j, :], q1[:, i]).reshape(
	                3,
	            )


	    for i in range(q1.shape[1]):
	        arrow_ = rotation(arrow,q1[:,i]).reshape(3,)
	        self.ax.plot(0,0,0,'ko')
	        self.ax.plot(framex_[:, 0, i], framex_[:, 1, i], framex_[:, 2, i], "r")
	        self.ax.plot(framey_[:, 0, i], framey_[:, 1, i], framey_[:, 2, i], "g")
	        self.ax.plot(framez_[:, 0, i], framez_[:, 1, i], framez_[:, 2, i], "b")
	        self.ax.plot(connect_zy_[:, 0, i], connect_zy_[:, 1, i], connect_zy_[:, 2, i], "k")
	        self.ax.plot(connect_zy1_[:, 0, i],connect_zy1_[:, 1, i], connect_zy1_[:, 2, i], "k")
	        self.ax.plot(self.connect_zI[:, 0], self.connect_zI[:, 1], self.connect_zI[:, 2], "k")
	        self.ax.plot(self.connect_zI1[:, 0], self.connect_zI1[:, 1], self.connect_zI1[:, 2], "k")
	        self.ax.plot(connect_xy_[:, 0, i], connect_xy_[:, 1, i], connect_xy_[:, 2, i], "k")
	        self.ax.plot(connect_xy1_[:, 0, i], connect_xy1_[:, 1, i], connect_xy1_[:, 2, i], "k")
	        a = Arrow3D([0, arrow_[0]], [0, arrow_[1]], [0, arrow_[2]], mutation_scale=10,
	            lw=1, arrowstyle="-|>", color="k")
	        self.ax.add_artist(a)
	        self.ax.set_xlabel("X")
	        self.ax.set_ylabel("Y")
	        self.ax.set_zlabel("Z")
	        self.ax.set_xlim([-2,2])
	        self.ax.set_ylim([-2,2])
	        self.ax.set_zlim([-2,2])
	        self.ax.plot(body_[:,0,i],body_[:,1,i],body_[:,2,i],'k')
	       
	    
	        plt.pause(0.001)  # pause avec duree en secondes
	        if i < q1.shape[1] - 1:
	        	a.remove()
	        	self.ax.clear()
	        plt.draw()
	    
	    self.q1_ = q1[:,-1]
	    self.q2_ = q2[:,-1]
	    self.q3_ = q3[:,-1]
	    
	            
	    return q1, q2, q3

	
def handle_click():
    phi1 = float(phi.get())
    theta1 = float(theta.get())
    psi1 = float(psi.get())
    time1 = float(time.get())
    Ts1 = float(Ts.get())
    
    d2r = np.pi/180
    a1,a2,a3 = g.animation(phi1*d2r,theta1*d2r,psi1*d2r,time1,Ts1)
    quat1.append(np.block([Ts1*np.ones((a1.shape[1],1)),a1.T]))
    quat2.append(np.block([Ts1*np.ones((a2.shape[1],1)),a2.T]))    
    quat3.append(np.block([Ts1*np.ones((a3.shape[1],1)),a3.T]))
    angles = quaternion2Euler(a1.T)
    Euler321rad.append(angles)
    
    print("---------------")
    print("Euler Angles 321 (deg): ", angles[-1,:]/d2r)
    print("Quaternion 0123: ", a1[:,-1])
    
    
    
if __name__ == "__main__":
    
    g = gimbals()
    
    window = tk.Tk()
    window.title("Gimbals")
    frame = tk.Frame(master=window)
    frame.grid(row=0, column=0, padx=0)
    
    w = 15 #width for the entries
    v1 = tk.StringVar(window, value='0')
    phi = tk.Entry(master=frame,width=w,textvariable=v1)
    phi.grid(row=1, column=0)
    lbl_phi = tk.Label(master=frame, text="Angle for the inner gimbal in (degrees)",)
    lbl_phi.grid(row=1, column=1)
    
    v2 = tk.StringVar(window, value='0')
    theta = tk.Entry(master = frame, width=w,textvariable=v2)
    theta.grid(row=2, column=0)
    lbl_theta = tk.Label(master=frame, text="Angle for the middle gimbal in (degrees)")
    lbl_theta.grid(row=2, column=1)
    
    v3 = tk.StringVar(window, value='0')
    psi = tk.Entry(master=frame, width=w,textvariable=v3)
    psi.grid(row=3, column=0)
    lbl_psi = tk.Label(master=frame, text="Angle for the outer gimbal in (degrees)")
    lbl_psi.grid(row=3, column=1)
    
    t = tk.StringVar(window, value='5') 
    time = tk.Entry(master=frame, width=w, textvariable=t)
    time.grid(row=4, column=0)
    lbl_time = tk.Label(master=frame, text="Simulation time (s)")
    lbl_time.grid(row=4, column=1)
    
    ts = tk.StringVar(window, value='0.1')
    Ts = tk.Entry(master=frame, width=w, textvariable=ts)
    Ts.grid(row=5, column=0)
    lbl_Ts = tk.Label(master=frame, text="Sample time (s)")
    lbl_Ts.grid(row=5, column=1)
    
    g.animation(0,0,0,0.1,0.1)
    
    
    button = tk.Button(master = frame, text="Run", command=handle_click, width=10)
    button.grid(row=7,column=0)
    window.mainloop()
    
    
	
	





