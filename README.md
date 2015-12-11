# ŌvSim

Tools for the simulation and visualization of dynamic follicle development within the mammalian ovary have not been available. We hypothesized that establishing a simple set of rules for i) follicle growth activation, ii) granulosa cell proliferation, and iii) follicle survival could provide the necessary starting points for a rudimentary simulation of stochastic follicle behavior over time. 

We produced a function in the R language, providing the ability to conduct simulations of ovarian follicle population development that vary based on the above user-specified parameters (or inputs). 


***JAY: A SHORT PARAGRAPH ABOUT R SPECIFICS. A PACKAGE? A SCRIPT? BRIEF INSTRUCTIONS ABOUT HOW TO USE GO HERE.***


To our surprise, the simple probability model can produce remarkably accurate representations of follicle population dynamics, closely matching the biologically observed number of surviving follicles (and thus an estimate of ovulated eggs) over time. Although this does not prove that the apparently complex process of follicle population dynamics is simple, the results show that a relatively simple probability-dependent process is consistent with and could help us better understand the process of follicle development in nature.

ŌvSim was designed using known biological parameters of ovarian follicles while allowing users to set certain variables. Once the script is activated, a numerical matrix is populated with values corresponding to the starting number of granulosa cells in individual simulated follicles. The matrix is initially populated with primordial follicles, randomly-generated such that they can contain one, two, or three pregranulosa cells.

For a truncated example, consider *n* matrix entries (containing 1, 2, or 3), representing *n* growth arrested ovarian "follicles" at start. Upon activating the simulation, a subset of follicles can growth activate at each step, while the remainder will remain growth-arrested and maintain an unchanged number of pregranulosa cells. Instead, growing follicles double in population at each step, minus a user-controlled factor of granulosa cell death. At each step, follicles thus grow exponentially (minus the cell death factor) or die. If a follicle that had begun to grow commits to death (atresia), its granulosa cell number is set to zero.

In ŌvSim, the starting number of follicles in the ovary (*NF*), the number of days of time (*ND*) to run the simulation, and the length of the ovulatory cycle (*cyclength*) can all be specified. We set the number of mouse ovarian follicles to 3000, including 2250 primordial follicles for most of our studies. Ovulatory cycle length for mice is set at 4, 4.5, or 5 days. The script then continues to loop with "daily" probability calculations and operations upon each follicle entry in the matrix. A flow chart of the operations upon each matrix entry is shown in Figure 1. 

![Fig1](OvSim_Fig1.jpg)

Variables are entered into the *follicle* function as shown in the following example. 

	follicle <- function(NF = 3000,
         ND = 420,
         IGP = 300,
         phold = 0.995,
         cond.pdub = 0.9,
         pcelllive = 0.8,
         cyclength = 5,
         ejectnum = 50000,
         puberty = TRUE,
         verbose = TRUE,
         pdfname = NA)
	{


Here, 3000 total follicles are present at the start, 2250 of which are primordial (1-3 granulosa cells), and 750 are small growing follicles randomly modeled to contain up to 20,000 granulosa cells. The estrus cycle length is 5 days. The probabilities shown next to the arrows in the flow chart in Figure 1 match these settings. Simulations are run for 420 days (14 months), approximately the fertile lifespan of C57Bl/6 mice fed *ad libitum*.

ŌvSim also generates .PDF output, the defaults of which are included as follows. 
Panel 2A shows the trend of decline of the primordial pool over time after execution of the simulation 1000 times. The plot represents the 1st and 99th percentiles of simulation data (A, gray hatched area). Individual data points for actual counts of C57Bl/6 mouse follicles in histological sections at 6 weeks, 8 months, and 12 months are overlaid with simulated data (circles). 

![2A](2A.png)

Panel 2B is a plot of the growth and death of individual follicles that die within the 420 days of simulated time. Granulosa cell number is represented by the dashed lines, and the time (and follicle ``size") of death is indicated by the letter ``D."

![2B](2B.png)

Last, Panel 2C is a histogram plot of the distribution of follicles that survive to ovulatory size, grouped in 4.5 day increments equivalent to the modeled estrus (e.g., ovulatory) cycle length. The number of eggs available for ovulation each cycle are therefore depicted. ŌvSim also outputs these data in Panel C in text format, useful for finer analyses of ovulatory follicle number. 

![2C](2C.png)

