# modified from: https://gist.githubusercontent.com/payoung/6087046/raw/a73cccefaed44d5899319469e30b6b325dbda665/tsp_plot.py
import os
import matplotlib.pyplot as plt
from pathlib import Path
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (18, 14),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
plt.rcParams.update(params)

# global vars
LOG_PATH = "logs/002.log"
CLEAR_LOG_FILE = False
INCLUDE_RED_ARROWS = True

def plotTSP(city_tours):
    for city_name, city, tours in city_tours:
        for idx, (tour_name, tour) in enumerate(tours):
            x = []; y = []
            for i in tour:
                x.append(city[i][0])
                y.append(city[i][1])
            plt.plot(x, y, 'co')
            a_scale = float(max(x))/float(300)

            if INCLUDE_RED_ARROWS:
                for i, (tour_name_i, tour_i) in enumerate(tours[idx+1:]):
                    num_iters = len(tours) - idx
                    xi = []; yi = [];
                    for j in tour_i:
                        xi.append(city[j][0])
                        yi.append(city[j][1])
                    plt.arrow(xi[-1], yi[-1], (xi[0] - xi[-1]), (yi[0] - yi[-1]),
                            head_width = a_scale, color = 'r',
                            length_includes_head = True, ls = 'dashed',
                            width = 0.001/float(num_iters))
                    for i in range(0, len(x) - 1):
                        plt.arrow(xi[i], yi[i], (xi[i+1] - xi[i]), (yi[i+1] - yi[i]),
                                head_width = a_scale, color = 'r', length_includes_head = True,
                                ls = 'dashed', width = 0.001/float(num_iters))

            # Draw the primary path for the TSP problem
            plt.arrow(x[-1], y[-1], (x[0] - x[-1]), (y[0] - y[-1]), head_width = a_scale,
                    color ='g', length_includes_head=True)
            for i in range(0,len(x)-1):
                plt.arrow(x[i], y[i], (x[i+1] - x[i]), (y[i+1] - y[i]), head_width = a_scale,
                        color = 'g', length_includes_head = True)

            #Set axis too slitghtly larger than the set of x and y
            plt.xlim(0, max(x)*1.1)
            plt.ylim(0, max(y)*1.1)

            filename =  f"{city_name} -- t:{tour_name}"
            plt.title(filename)
            Path(os.path.join(LOG_PATH+"_pics", city_name)).mkdir(parents=True, exist_ok=True)
            plt.savefig(os.path.join(os.path.join(LOG_PATH+"_pics", city_name), filename + ".png"), dpi=300)
            # plt.show()
            plt.close()

        # for tour_name, tour in tours:
        #     plt.plot(x, y, 'co')

        #     # Set a scale for the arrow heads (there should be a reasonable default for this??)
        #     a_scale = float(max(x))/float(300)

        #     # for i, tour_iteration in enumerate(tour[1:]):
        #     #     xi = []; yi = [];
        #     #     for j in tour_iteration:
        #     #         xi.append(points[j][0])
        #     #         yi.append(points[j][1])

        #     #     plt.arrow(xi[-1], yi[-1], (xi[0] - xi[-1]), (yi[0] - yi[-1]),
        #     #             head_width = a_scale, color = 'r',
        #     #             length_includes_head = True, ls = 'dashed',
        #     #             width = 0.001/float(num_iters))
        #     #     for i in range(0, len(x) - 1):
        #     #         plt.arrow(xi[i], yi[i], (xi[i+1] - xi[i]), (yi[i+1] - yi[i]),
        #     #                 head_width = a_scale, color = 'r', length_includes_head = True,
        #     #                 ls = 'dashed', width = 0.001/float(num_iters))

        # # Draw the primary path for the TSP problem
        # plt.arrow(x[-1], y[-1], (x[0] - x[-1]), (y[0] - y[-1]), head_width = a_scale,
        #         color ='g', length_includes_head=True)
        # for i in range(0,len(x)-1):
        #     plt.arrow(x[i], y[i], (x[i+1] - x[i]), (y[i+1] - y[i]), head_width = a_scale,
        #             color = 'g', length_includes_head = True)

        # #Set axis too slitghtly larger than the set of x and y
        # plt.xlim(0, max(x)*1.1)
        # plt.ylim(0, max(y)*1.1)

        # filename =  f"{city_name} -- t:{tour_name}"
        # plt.title(filename)
        # plt.savefig(os.path.join(LOG_PATH+"_pics", filename + ".png"), dpi=300)
        # # plt.show()
        # plt.close()


if __name__ == '__main__':
    city_tours = [] # [ ..., [city_name, city, [  [tour1_name, tour1],  [tour2_name, tour2]]  ], ...]
    Path(LOG_PATH+"_pics").mkdir(parents=True, exist_ok=True)
    with open(LOG_PATH, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("cities_name="):
                city_name = line.replace("cities_name=", "")
                city_tours.append([city_name])
            if line.startswith("cities="):
                city = eval(line.replace("cities=", ""))
                city_tours[-1].append(city)
                city_tours[-1].append([])
            if line.startswith("tour_name="):
                tour_name = line.replace("tour_name=", "")
                city_tours[-1][-1].append([tour_name])
            if line.startswith("tour="):
                tour = eval(line.replace("tour=", ""))
                city_tours[-1][-1][-1].append(tour)
    if CLEAR_LOG_FILE:
        open(LOG_PATH, "w").write("")
    # print(city_tours)

    plotTSP(city_tours)

