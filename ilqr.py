import rospy
import numpy as np 
import copy
import action, dataSample
from sklearn.linear_model import LinearRegression
from autograd import grad, jacobian


class iLQG:
	def __init__(self, umax, state_dim, action_dim, pred_time = 50, n = 5):


		# take sequential sample from real robot
		self.x_seq = [np.zeros(self.state_dim) for _ in range(pred_time+1)]
		self.x_seq_mean = copy.deepcopy(self.x_seq)
		self.u_seq = [np.zeros(self.action_dim) for _ in range(pred_time)]

		#self.cnt = 0 # callback counter
		#argument n: each policy repete times
		self.n = n
		self.pred_time = pred_time
		self.action_dim = action_dim
		self.state_dim = state_dim
		self.umax = umax
		self.v = [0.0 for _ in range(pred_time+1)]
		self.v_x = [np.zeros(state_dim) for _ in range(pred_time+1)]
		self.v_xx = [np.zeros((state_dim, state_dim)) for _ in range(pred_time+1)]
		self.lamb = 0.99  #dual variable

		self.action_mu = [np.zeros(self.action_dim) for _ in range(pred_time)]
		self.action_sigma = [ np.ones((action_dim, action_dim)) for _ in range(pred_time)]
		self.action_mu_prev = [np.zeros(self.action_dim) for _ in range(pred_time)]
		self.action_sigma_prev = [ np.ones((action_dim, action_dim)) for _ in range(pred_time)]

		self.f_x = [np.zeros((state_dim, state_dim)) for _ in range(pred_time)]
		self.f_u = [np.zeros((state_dim, action_dim)) for _ in range(pred_time)]
		self.lf = lambda x: 0.5*np.sum(np.square(x)) # lf(x) 
		self.lf_x = grad(self.lf)
		self.lf_xx = jacobian(self.lf_x)
		# l_prime(x,u) = (1/lambda)*l(x,u) - log(p(u)) 

		self.kk = [np.zeros((self.action_dim, state_dim)) for _ in range(pred_time)]
		self.k = [np.zeros(self.action_dim) for _ in range(pred_time)]


	def dynamic_fit(self, x_, x, u):
		LR = LinearRegression(fit_intercept=False)
		A = []
		B = []
		for t in range(self.pred_time):

			#x[t] and u[t] is both np.array type, shape: num_sample*(state_dim + action_dim)
			input_x = np.hstack((x[t], u[t]))	
			
			# x_[t] is np.array type. Shape: num_sample * state_dim
			reg = LR.fit(input_x, x_[t])
			A.append( reg.coef_[:, :self.state_dim])
			B.append( reg.coef_[:, self.state_dim:])
		self.f_x = copy.deepcopy(A)
		self.f_u = copy.deepcopy(B)


	def eta_update(self, x_seq, u_seq):
		D_kl=0
		mean_cross_entropy = 0
		step_size = 0.001

		for t in range(self.pred_time):
			for i in range(u_seq[t].shape[0]):				
				u = u_seq[t,i,:]
				mean_cross_entropy += (-self.action_dim/2)*np.log(2*np.pi)-0.5*np.log(np.linalg.det(self.action_sigma)) - \
				 0.5*np.matmul(np.matmul((u - self.action_mu).T, np.linalg.inv(self.action_sigma)), (u-self.action_mu)) +\
				 (self.action_dim/2)*np.log(2*np.pi) + 0.5*np.log(np.linalg.det(self.action_sigma_prev)) + \
				 0.5*np.matmul(np.matmul((u - self.action_mu_prev).T, np.linalg.inv(self.action_sigma_prev)), (u-self.action_mu_prev))
			D_kl += mean_cross_entropy/u_seq[t].shape[0]
			mean_cross_entropy = 0

		self.lamb = self.lamb + step_size*D_kl
		

	def backward(self, x, u):

		self.running_cost = lambda x, u: (1/self.lamb)*np.sum(np.square(x)) + \
		 (self.action_dim/2)*np.log(2*np.pi) + 0.5*np.log(np.linalg.det(self.action_sigma)) + \
		 0.5*np.matmul(np.matmul((u - self.action_mu).T, np.linalg.inv(self.action_sigma)), (u-self.action_mu))

		self.l_x = grad(self.running_cost, 0)
		self.l_u = grad(self.running_cost, 1)
		self.l_xx = jacobian(self.l_x, 0)
		self.l_uu = jacobian(self.l_u, 1)
		self.l_ux = jacobian(self.l_u, 0)

		self.v[-1] = self.lf(x[-1])
		self.v_x[-1] = self.lf_x(x[-1])
		self.v_xx[-1] = self.lf_xx(x[-1])
		k_seq = []
		kk_seq = []
		for t in range(self.pred_time - 1, -1, -1):
			q_x = self.l_x(x[t], u[t]) + np.matmul(self.f_x[t].T, self.v_x[t + 1])
			q_u = self.l_u(x[t], u[t]) + np.matmul(self.f_u[t].T, self.v_x[t + 1])
			q_xx = self.l_xx(x[t], u[t]) + np.matmul(np.matmul(self.f_x[t].T, self.v_xx[t+1]), self.f_x[t])
			tmp = np.matmul(self.f_u[t].T, self.v_xx[t+1])
			q_uu = self.l_uu(x[t], u[t]) + np.matmul(tmp, self.f_u[t])
			q_ux = self.l_ux(x[t], u[t]) + np.matmul(tmp, self.f_x[t])
			inv_q_uu = np.linalg.inv(q_uu)
			self.action_sigma = inv_q_uu
			k = -np.matmul(inv_q_uu, q_u)
			kk = -np.matmul(inv_q_uu, q_ux)
			dv = 0.5*np.matmul(q_u, k)
			self.v[t] += dv
			self.v_x[t] = q_x - np.matmul(np.matmul(q_u, inv_q_uu), q_ux)
			self.v_xx[t] = q_xx + np.matmul(q_ux.T, kk)
			k_seq.append(k)
			kk_seq.append(kk)
		k_seq.reverse()
		kk_seq.reverse()
		self.kk = kk_seq
		self.k = k_seq
		return kk_seq, k_seq


	def actionPub(self):

		actionValue = self.policy_to_msg()
		rospy.loginfo(actionValue)
		self.pub.publish(actionValue)
		self.rate.sleep()


	def callback(self, data):
		cnt = 0
		x_seq = []
		x = []
		x__seq = []
		u_seq = []
		u = []

		for d in data.data:
			if cnt < self.state_dim:
				x.append(d)
				cnt += 1
				if cnt >= self.state_dim:
					x_seq.append(np.array(x))
					x = []
			elif self.state_dim <= cnt and cnt < self.state_dim+self.action_dim:
				u.append(d)
				cnt += 1
				if cnt >= self.state_dim+self.action_dim:
					u_seq.append(np.array(u))
					u = []
					cnt = 0

# This is only one episode !!! should have 5 episodes samples????
		if self.x_seq[0] == np.zeros(self.state_dim) or self.x_seq[0].shape[0] == self.n:
			self.x_seq = copy.deepcopy(x_seq)
			self.u_seq = copy.deepcopy(u_seq)

		else:
			for i in range(len(x_seq)):
				self.x_seq[i] = np.vstack((self.x_seq[i],x_seq[i]))

			for j in range(len(u_seq)):
				self.u_seq[j] = np.vstack((self.u_seq[j],u_seq[j]))
		
		if self.x_seq[0].shape[0] == self.n and self.u_seq[0].shape[0] == self.n:

			self.action_mu_prev = self.action_mu
			self.action_sigma_prev = self.action_sigma
			self.action_mu = []
			self.action_sigma = []	

			for k in self.u_seq:
				self.action_mu.append(np.mean(k, axis=0))
				self.action_sigma.append(np.cov(k.T))

			self.x_seq_mean = []
			for l in self.x_seq:
				self.x_seq_mean.append(np.mean(l, axis=0))

	
	def run_iLQG(self):

		while not rospy.is_shutdown():

			for i in range(self.n): #policy run n times with noise
				self.actionPub()
				#rospy.Duration(5).sleep() # Duration 5 seconds to take sample trajectory
				rospy.sleep(5)
			self.eta_update(self.x_seq[:], self.u_seq[:])
			self.dynamic_fit(self.x_seq[1:], self.x_seq[:-1], self.u_seq[:])
			self.backward(self.x_seq_mean[:], self.action_mu[:])

		#if not rospy.is_shutdown():
			#rospy.spin()


	def policy_to_msg(self):
		msg = action()
		msg.dX = self.state_dim
		msg.dU = self.action_dim
		msg.K_t = []
		msg.k_t = []
		msg.x_seq_prev =[]
		msg.u_seq_prev = []
		for kk in self.kk:
			msg.K_t.extend(kk.reshape(self.state_dim*self.action_dim).tolist())

		for i in range(len(self.k)):
			k = self.k[i] + np.random.multivariate_normal([0,0,0,0,0], self.action_sigma[i])
			msg.k_t.extend(k.reshape(self.action_dim).tolist())

		for ac in self.action_mu:
			msg.u_seq_prev.extend(ac.reshape(self.action_dim).tolist())

		for x in self.x_seq_mean:
			msg.x_seq_prev.extend(x.reshape(self.state_dim).tolist())

		return msg 





