import numpy as np
import copy
from sklearn.linear_model import LinearRegression
from autograd import grad, jacobian


class iLQG:
	def __init__(self, umax, state_dim, action_dim, pred_time = 50, n = 5):
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

		# Sample sequence
		self.x_seq = []
		self.x_seq_mean = []
		self.x_seq_prev = []
		self.u_seq = []

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
			A.append(reg.coef_[:, :self.state_dim])
			B.append(reg.coef_[:, self.state_dim:])
		self.f_x = copy.deepcopy(A)
		self.f_u = copy.deepcopy(B)


	def eta_update(self, x_seq, u_seq):
		D_kl = 0.0
		step_size = 0.001
		eps = 1.0

		for t in range(self.pred_time):
			inv_action_sigma_prev = np.linalg.inv(self.action_sigma_prev[t])
			action_diff = self.action_mu_prev[t] - self.action_mu[t]

			D_kl += 0.5*(np.log(np.linalg.det(self.action_sigma_prev[t]) / np.linalg.det(self.action_sigma[t])) -\
			             self.action_dim + np.dot(inv_action_sigma_prev, self.action_sigma[t]).trace() + \
			             np.dot(np.dot(action_diff.T, inv_action_sigma_prev), action_diff))

		self.lamb = self.lamb + step_size * (D_kl - eps)
		print("The total Kl divergence now is updated to : {}".format(D_kl))
		print("lamb has been updated to {}".format(self.lamb))
		

	def backward(self, x, u):

		k_seq = []
		kk_seq = []
		weight = 1.0

		self.v[-1] = self.lf(x[-1])
		self.v_x[-1] = self.lf_x(x[-1])
		self.v_xx[-1] = self.lf_xx(x[-1])

		for t in range(self.pred_time - 1, -1, -1):
			# print(" u[t]'s type is :{}".format(type(u[t])))
			self.running_cost = lambda x_t, u_t: (1 / self.lamb) * np.sum(np.square(x_t)) + \
			                                     weight * (self.action_dim / 2) * np.log(2 * np.pi) + weight * 0.5 * np.log(np.linalg.det(self.action_sigma_prev[t])) + \
			                                     weight * 0.5 * np.dot(np.dot((u_t - self.action_mu_prev[t]).T, np.linalg.pinv(self.action_sigma_prev[t])),(u_t - self.action_mu_prev[t]))
			# self.running_cost = lambda x_t, u_t: (1/self.lamb)*np.sum(np.square(x_t)) +(1/self.lamb)*np.sum(np.square(u_t))
			# self.running_cost = lambda x_t, u_t: (1/self.lamb)*np.sum(np.square(x_t)) + 0.5*np.dot(np.dot((u_t - self.action_mu_prev[t]).T, np.linalg.pinv(self.action_sigma_prev[t])), (u_t-self.action_mu_prev[t]))
			self.l_x = grad(self.running_cost, 0)
			self.l_u = grad(self.running_cost, 1)
			self.l_xx = jacobian(self.l_x, 0)
			self.l_uu = jacobian(self.l_u, 1)
			self.l_ux = jacobian(self.l_u, 0)

			q_x = self.l_x(x[t], u[t]) + np.matmul(self.f_x[t].T, self.v_x[t + 1])
			q_u = self.l_u(x[t], u[t]) + np.matmul(self.f_u[t].T, self.v_x[t + 1])
			q_xx = self.l_xx(x[t], u[t]) + np.matmul(np.matmul(self.f_x[t].T, self.v_xx[t + 1]), self.f_x[t])
			tmp = np.matmul(self.f_u[t].T, self.v_xx[t + 1])
			q_uu = self.l_uu(x[t], u[t]) + np.matmul(tmp, self.f_u[t])
			q_ux = self.l_ux(x[t], u[t]) + np.matmul(tmp, self.f_x[t])

			self.action_sigma_prev[t] = self.action_sigma[t]
			inv_q_uu = np.linalg.pinv(q_uu)

			self.action_sigma[t] = inv_q_uu
			k = -np.matmul(inv_q_uu, q_u)
			kk = -np.matmul(inv_q_uu, q_ux)
			dv = 0.5 * np.matmul(q_u, k)
			self.v[t] += dv
			self.v_x[t] = q_x - np.matmul(np.matmul(q_u, inv_q_uu), q_ux)
			self.v_xx[t] = q_xx + np.matmul(q_ux.T, kk)
			k_seq.append(k)
			kk_seq.append(kk)

		k_seq.reverse()
		kk_seq.reverse()
		self.kk = kk_seq
		self.k = k_seq
		# print('backward action sigma:{}'.format(self.action_sigma))
		# print('backward action K:{}'.format(self.k))

		return kk_seq, k_seq


	def storeSample(self, x_cur, u_cur):
		assert x_cur[0].shape[0] == self.state_dim
		assert u_cur[0].shape[0] == self.action_dim

		if len(self.x_seq) == 0:
			for x in x_cur:
				self.x_seq.append(x)

			for u in u_cur:
				self.u_seq.append(u)
		else:
			for i in range(len(x_cur)):
				self.x_seq[i] = np.vstack((self.x_seq[i], x_cur[i]))

			for i in range(len(u_cur)):
				self.u_seq[i] = np.vstack((self.u_seq[i], u_cur[i]))

			if self.x_seq[0].shape[0] == self.n and len(self.x_seq[0].shape) == 2:
				self.x_seq_prev = copy.deepcopy(self.x_seq_mean)
				for i in range(len(self.x_seq)):
					self.x_seq_mean[i] = np.mean(self.x_seq[i], axis = 0)

				self.action_mu_prev = copy.deepcopy(self.action_mu)
				for i in range(len(self.u_seq)):
					self.action_mu[i] = np.mean(self.u_seq[i], axis = 0)


	def totalCost(self):
		immediateCost = 0
		weight = 1

		for t in range(self.pred_time):
			immediateCost += (1/self.lamb)*np.sum(np.square(self.x_seq_mean[t])) +\
			                 weight * (self.action_dim/2)*np.log(2*np.pi) + weight * 0.5*np.log(np.linalg.det(self.action_sigma_prev[t])) + \
			                 weight * 0.5*np.dot(np.dot((self.action_mu[t] - self.action_mu_prev[t]).T, np.linalg.pinv(self.action_sigma_prev[t])), (self.action_mu[t]-self.action_mu_prev[t]))
		return immediateCost


	def clearSample(self):
		self.x_seq = []
		self.u_seq = []


	def update(self):
		self.eta_update(self.x_seq[:], self.u_seq[:])
		self.dynamic_fit(self.x_seq[1:], self.x_seq[:-1], self.u_seq[:])
		self.backward(self.x_seq_mean[:], self.action_mu[:])


	def policy(self, state, index):
		if len(self.x_seq_mean) == 0:
			return np.random.multivariate_normal(np.zeros(self.action_dim), self.action_sigma[index])
		else:
			torque = self.kk[index] * (state - self.x_seq_mean[index]) + self.k[index] + self.action_mu[index] + \
			         np.random.multivariate_normal(np.zeros(self.action_dim), self.action_sigma[index])

		return torque








