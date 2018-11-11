import ilqr
import dynamic




if __name__ == '__main__':

	agent = ilqr.iLQG(umax = 100, state_dim = 10, action_dim = 5, pred_time = 50, n = 5 )
	agent.run_iLQG()
	

		


