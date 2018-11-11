import ilqr
import dynamic
import numpy as np
import copy



if __name__ == '__main__':

	epoch = 10000 # Trainning time
	_pred_time = 50 # The length of each episode
	_n = 5  # The number of running each policy

	agent = ilqr.iLQG(umax = 10000, state_dim = 10, action_dim = 5, pred_time = _pred_time, n = _n )
	env = dynamic.dynamic(dof=5, delta_t=0.02)


	for epo in range(epoch):

		for n in range(_n):

			x_sequence = []
			u_sequence = []
			env.init_sys()  # Initialize the system
			assert  env.pos[0] == 0

			for pred in range(_pred_time):
				x_cur = np.hstack((env.pos, env.vel))
				x_sequence.append(x_cur)
				torque = agent.policy(x_cur, pred)
				u_sequence.append(torque)
				env.transition(torque)

			x_sequence.append(np.hstack((env.pos, env.vel)))
			agent.storeSample(x_sequence, u_sequence)

		## Start learning
		agent.update()
		agent.clearSample()


		


