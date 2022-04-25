# Simulation of a double pendulum with 2 equal masses

# inspired by https://matplotlib.org/stable/gallery/animation/double_pendulum.html#sphx-glr-gallery-animation-double-pendulum-py


import matplotlib.pyplot as plt
import matplotlib.animation
import scipy.integrate
import numpy as np


class DoublePendulum:
    '''
    define the problem: double pendulum
    all parameters are defined
    '''

    def __init__(self, parameters_dict):
        self.dimension_rhs = 4
        self.initial_state = np.zeros(self.dimension_rhs)
        for i in range(1, int(self.dimension_rhs / 2 + 1)):
            self.initial_state[2 * i - 2] = parameters_dict["initial_angle%d" % i] * np.pi / 180
            self.initial_state[2 * i - 1] = parameters_dict["initial_velocity%d" % i] * np.pi / 180
        self.time = 0.0
        self.mass = parameters_dict["mass"]
        self.length = parameters_dict["length"]
        self.gravity = parameters_dict["gravity"]
        self.step_size_animation = parameters_dict["time_step_size_animation"]

    def rhs(self, time, state):
        z = np.zeros_like(self.initial_state)
        z[0] = state[1]
        z[2] = state[3]
        delta = state[2] - state[0]
        M = np.array([[2, np.cos(delta)], [np.cos(delta), 1]])
        v = np.array([-2 * self.gravity / self.length * np.sin(state[0]) + np.sin(delta) * state[3] ** 2,
                      -self.gravity / self.length * np.sin(state[2]) - np.sin(delta) * state[1] ** 2])
        [z[1], z[3]] = np.linalg.solve(M, v)
        return z


def animate_pendulum(problem, cartesian_coord, dt):
    # creates an animation of the pendulum
    # input: solution (position vector) in cartesian coordinates [x1,y1,x2,y2]
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(autoscale_on=True, xlim=(-2 * problem.length, 2 * problem.length),
                         ylim=(-2.1 * problem.length, 1.))
    ax.set_aspect('equal')
    ax.grid()

    # creating a line plot with circles at the end
    line, = ax.plot([], [], 'o-', lw=2)
    trace, = ax.plot([], [], '-', lw=1, ms=2)  # plotting the pendulum's trace

    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
    history_x, history_y = [], []

    def update(i):
        # update the position of the pendulums
        x1 = cartesian_coord[0][0]  # pendulum 1
        y1 = cartesian_coord[0][1]
        x2 = cartesian_coord[1][0]  # pendulum 2
        y2 = cartesian_coord[1][1]
        thisx = [0, x1[i], x2[i]]
        thisy = [0, y1[i], y2[i]]
        if i == 0:
            history_x.clear()
            history_y.clear()

        history_x.append(thisx[2])
        history_y.append(thisy[2])

        line.set_data(thisx, thisy)
        trace.set_data(history_x, history_y)
        time_text.set_text(time_template % (i * dt))
        return line, trace, time_text

    ani = matplotlib.animation.FuncAnimation(
        fig, update, len(cartesian_coord[0][0]), interval=dt * 1000, blit=True)
    f = r"animation.gif"
    writergif = matplotlib.animation.PillowWriter(fps=60)
    ani.save(f, writer=writergif)
    plt.show()


def main():
    # import parameters from file txt
    experiment = int(input("Insert test number: "))
    test_file = open('test_files/doublePendulum_test%d.txt' % experiment, 'r')
    header_line = next(test_file)
    initial_state, parameters_dict, problem, t_span = read_test_file(test_file)

    solution = scipy.integrate.solve_ivp(fun=problem.rhs, t_span=t_span, y0=initial_state, first_step=1e-04, atol=1e-6,
                                         rtol=1e-6, dense_output=True)
    dt = float(parameters_dict["time_step_size_animation"])
    t_dense = np.arange(t_span[0], t_span[1], dt)
    dense_solution = solution.sol(t_dense)

    # plot time-solution-diagram
    if int(parameters_dict["time_solution_diagram"]) == 1:
        plot_solution_over_time(dense_solution, t_dense)

    # x-position-velocity diagrams
    if int(parameters_dict["position_velocity_diagrams"]) == 1:
        pos = position_velocity_diagrams(problem, dense_solution)

    # animate pendulum
    cart_coord = polar_to_cartesian(problem, dense_solution)
    print(energy_test(problem, dense_solution))

    animation = animate_pendulum(problem, cart_coord, dt)


def plot_solution_over_time(dense_solution, t_dense):
    plt.figure()
    ax = plt.gca()
    for i in range(int(len(dense_solution) / 2)):
        ax.plot(t_dense, dense_solution[2 * i], label="angle%d" % (i + 1))
    for i in range(int(len(dense_solution) / 2)):
        ax.plot(t_dense, dense_solution[2 * i + 1], label="velocity%d" % (i + 1))
    plt.xlabel('Time'), plt.ylabel('Angle and angular velocity')
    ax.legend()


def position_velocity_diagrams(problem, solution_vector):
    # x-position-velocity diagrams
    plt.figure()
    pos = polar_to_cartesian(problem, solution_vector)
    dx1dt = problem.length * solution_vector[1] * np.cos(solution_vector[0])
    dx2dt = dx1dt + problem.length * solution_vector[3] * np.cos(solution_vector[2])
    plt.subplot(1, 2, 1), plt.plot(pos[0][0], dx1dt), \
    plt.xlabel(r'$x_1$'), plt.ylabel(r'$\dot{x}_1$')
    plt.subplot(1, 2, 2), plt.plot(pos[1][0], dx2dt), \
    plt.xlabel(r'$x_2$'), plt.ylabel(r'$\dot{x}_2$')
    plt.suptitle('Phase Plots, cartesian coordinates', fontsize=20)
    return pos


def polar_to_cartesian(problem, solution):
    x1 = problem.length * np.sin(solution[0])
    y1 = -problem.length * np.cos(solution[0])
    x2 = x1 + problem.length * np.sin(solution[2])
    y2 = y1 - problem.length * np.cos(solution[2])
    return [(x1, y1), (x2, y2)]


def polar_to_cartesian_velocity(problem, solution_vector):
    dx1dt = problem.length * solution_vector[1] * np.cos(solution_vector[0])
    dx2dt = dx1dt + problem.length * solution_vector[3] * np.cos(solution_vector[2])
    dy1dt = problem.length * solution_vector[1] * np.sin(solution_vector[0])
    dy2dt = dy1dt + problem.length * solution_vector[3] * np.sin(solution_vector[2])
    return [(dx1dt, dy1dt), (dx2dt, dy2dt)]


def read_test_file(test_file):
    parameters_dict = {}
    for line in test_file:
        names, values = line.split(' = ')
        parameters_dict[names] = float(values)
    problem = DoublePendulum(parameters_dict)
    initial_state = problem.initial_state  # [angle1,angular_velocity1,angle2,angular_velocity2]
    t_span = [0, int(parameters_dict["time_interval"])]  # =[t0,t_end] integration interval
    test_file.close()
    return initial_state, parameters_dict, problem, t_span


def energy_test(problem, solution_vector):
    # calculates the total energy of the system
    # solution_vector= [angle1,velocity1,angle2,velocity2]

    cartesian_solution = polar_to_cartesian(problem, solution_vector)
    y = np.array([cartesian_solution[i][1] for i in range(2)])
    cartesian_velocity = polar_to_cartesian_velocity(problem, solution_vector)
    dxdt = np.array([cartesian_velocity[i][0] for i in range(2)])
    dydt = np.array([cartesian_velocity[i][1] for i in range(2)])
    potential_energy = -problem.gravity * problem.mass * (sum(y))
    kinetic_energy = problem.mass / 2 * sum(dxdt ** 2 + dydt ** 2)
    total_energy = potential_energy + kinetic_energy
    return total_energy


if __name__ == '__main__':
    main()
