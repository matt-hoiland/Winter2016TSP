using System;
using System.Collections;
using System.Collections.Generic;
using System.Text;
using System.Drawing;
using System.Diagnostics;
using System.Timers;


namespace TSP
{

    class ProblemAndSolver
    {

        internal class TSPSolution
        {
            /// <summary>
            /// we use the representation [cityB,cityA,cityC] 
            /// to mean that cityB is the first city in the solution, cityA is the second, cityC is the third 
            /// and the edge from cityC to cityB is the final edge in the path.  
            /// You are, of course, free to use a different representation if it would be more convenient or efficient 
            /// for your data structure(s) and search algorithm. 
            /// </summary>
            public ArrayList
                Route;

            /// <summary>
            /// constructor
            /// </summary>
            /// <param name="iroute">a (hopefully) valid tour</param>
            public TSPSolution(ArrayList iroute)
            {
                Route = new ArrayList(iroute);
            }

            /// <summary>
            /// Compute the cost of the current route.  
            /// Note: This does not check that the route is complete.
            /// It assumes that the route passes from the last city back to the first city. 
            /// </summary>
            /// <returns></returns>
            public double costOfRoute()
            {
                // go through each edge in the route and add up the cost. 
                int x;
                City here;
                double cost = 0D;

                for (x = 0; x < Route.Count - 1; x++)
                {
                    here = Route[x] as City;
                    cost += here.costToGetTo(Route[x + 1] as City);
                }

                // go from the last city to the first. 
                here = Route[Route.Count - 1] as City;
                cost += here.costToGetTo(Route[0] as City);
                return cost;
            }
        }

        #region Private members 

        /// <summary>
        /// Default number of cities (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Problem Size text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int DEFAULT_SIZE = 25;

        /// <summary>
        /// Default time limit (unused -- to set defaults, change the values in the GUI form)
        /// </summary>
        // (This is no longer used -- to set default values, edit the form directly.  Open Form1.cs,
        // click on the Time text box, go to the Properties window (lower right corner), 
        // and change the "Text" value.)
        private const int TIME_LIMIT = 60;        //in seconds

        private const int CITY_ICON_SIZE = 5;


        // For normal and hard modes:
        // hard mode only
        private const double FRACTION_OF_PATHS_TO_REMOVE = 0.20;

        /// <summary>
        /// the cities in the current problem.
        /// </summary>
        private City[] Cities;
        /// <summary>
        /// a route through the current problem, useful as a temporary variable. 
        /// </summary>
        private ArrayList Route;
        /// <summary>
        /// best solution so far. 
        /// </summary>
        private TSPSolution bssf; 
        public TSPSolution BSSF
        {
            get { return bssf; } set { bssf = value; }
        }

        /// <summary>
        /// how to color various things. 
        /// </summary>
        private Brush cityBrushStartStyle;
        private Brush cityBrushStyle;
        private Pen routePenStyle;


        /// <summary>
        /// keep track of the seed value so that the same sequence of problems can be 
        /// regenerated next time the generator is run. 
        /// </summary>
        private int _seed;
        /// <summary>
        /// number of cities to include in a problem. 
        /// </summary>
        private int _size;

        /// <summary>
        /// Difficulty level
        /// </summary>
        private HardMode.Modes _mode;

        /// <summary>
        /// random number generator. 
        /// </summary>
        private Random rnd;

        /// <summary>
        /// time limit in milliseconds for state space search
        /// can be used by any solver method to truncate the search and return the BSSF
        /// </summary>
        private int time_limit;
        #endregion

        #region Public members

        /// <summary>
        /// These three constants are used for convenience/clarity in populating and accessing the results array that is passed back to the calling Form
        /// </summary>
        public const int COST = 0;           
        public const int TIME = 1;
        public const int COUNT = 2;
        
        public int Size
        {
            get { return _size; }
        }

        public int Seed
        {
            get { return _seed; }
        }
        #endregion

        #region Constructors
        public ProblemAndSolver()
        {
            this._seed = 1; 
            rnd = new Random(1);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed)
        {
            this._seed = seed;
            rnd = new Random(seed);
            this._size = DEFAULT_SIZE;
            this.time_limit = TIME_LIMIT * 1000;                  // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }

        public ProblemAndSolver(int seed, int size)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = TIME_LIMIT * 1000;                        // TIME_LIMIT is in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        public ProblemAndSolver(int seed, int size, int time)
        {
            this._seed = seed;
            this._size = size;
            rnd = new Random(seed);
            this.time_limit = time*1000;                        // time is entered in the GUI in seconds, but timer wants it in milliseconds

            this.resetData();
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Reset the problem instance.
        /// </summary>
        private void resetData()
        {

            Cities = new City[_size];
            Route = new ArrayList(_size);
            bssf = null;

            if (_mode == HardMode.Modes.Easy)
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble());
            }
            else // Medium and hard
            {
                for (int i = 0; i < _size; i++)
                    Cities[i] = new City(rnd.NextDouble(), rnd.NextDouble(), rnd.NextDouble() * City.MAX_ELEVATION);
            }

            HardMode mm = new HardMode(this._mode, this.rnd, Cities);
            if (_mode == HardMode.Modes.Hard)
            {
                int edgesToRemove = (int)(_size * FRACTION_OF_PATHS_TO_REMOVE);
                mm.removePaths(edgesToRemove);
            }
            City.setModeManager(mm);

            cityBrushStyle = new SolidBrush(Color.Black);
            cityBrushStartStyle = new SolidBrush(Color.Red);
            routePenStyle = new Pen(Color.Blue,1);
            routePenStyle.DashStyle = System.Drawing.Drawing2D.DashStyle.Solid;
        }

        #endregion

        #region Public Methods

        /// <summary>
        /// make a new problem with the given size.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode)
        {
            this._size = size;
            this._mode = mode;
            resetData();
        }

        /// <summary>
        /// make a new problem with the given size, now including timelimit paremeter that was added to form.
        /// </summary>
        /// <param name="size">number of cities</param>
        public void GenerateProblem(int size, HardMode.Modes mode, int timelimit)
        {
            this._size = size;
            this._mode = mode;
            this.time_limit = timelimit*1000;                                   //convert seconds to milliseconds
            resetData();
        }

        /// <summary>
        /// return a copy of the cities in this problem. 
        /// </summary>
        /// <returns>array of cities</returns>
        public City[] GetCities()
        {
            City[] retCities = new City[Cities.Length];
            Array.Copy(Cities, retCities, Cities.Length);
            return retCities;
        }

        /// <summary>
        /// draw the cities in the problem.  if the bssf member is defined, then
        /// draw that too. 
        /// </summary>
        /// <param name="g">where to draw the stuff</param>
        public void Draw(Graphics g)
        {
            float width  = g.VisibleClipBounds.Width-45F;
            float height = g.VisibleClipBounds.Height-45F;
            Font labelFont = new Font("Arial", 10);

            // Draw lines
            if (bssf != null)
            {
                // make a list of points. 
                Point[] ps = new Point[bssf.Route.Count];
                int index = 0;
                foreach (City c in bssf.Route)
                {
                    if (index < bssf.Route.Count -1)
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[index+1]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    else 
                        g.DrawString(" " + index +"("+c.costToGetTo(bssf.Route[0]as City)+")", labelFont, cityBrushStartStyle, new PointF((float)c.X * width + 3F, (float)c.Y * height));
                    ps[index++] = new Point((int)(c.X * width) + CITY_ICON_SIZE / 2, (int)(c.Y * height) + CITY_ICON_SIZE / 2);
                }

                if (ps.Length > 0)
                {
                    g.DrawLines(routePenStyle, ps);
                    g.FillEllipse(cityBrushStartStyle, (float)Cities[0].X * width - 1, (float)Cities[0].Y * height - 1, CITY_ICON_SIZE + 2, CITY_ICON_SIZE + 2);
                }

                // draw the last line. 
                g.DrawLine(routePenStyle, ps[0], ps[ps.Length - 1]);
            }

            // Draw city dots
            foreach (City c in Cities)
            {
                g.FillEllipse(cityBrushStyle, (float)c.X * width, (float)c.Y * height, CITY_ICON_SIZE, CITY_ICON_SIZE);
            }

        }

        /// <summary>
        ///  return the cost of the best solution so far. 
        /// </summary>
        /// <returns></returns>
        public double costOfBssf ()
        {
            if (bssf != null)
                return (bssf.costOfRoute());
            else
                return -1D; 
        }

        /// <summary>
        /// This is the entry point for the default solver
        /// which just finds a valid random tour 
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] defaultSolveProblem()
        {
            int i, swap, temp, count=0;
            string[] results = new string[3];
            int[] perm = new int[Cities.Length];
            Route = new ArrayList();
            Random rnd = new Random();
            Stopwatch timer = new Stopwatch();

            timer.Start();

            do
            {
                for (i = 0; i < perm.Length; i++)                                 // create a random permutation template
                    perm[i] = i;
                for (i = 0; i < perm.Length; i++)
                {
                    swap = i;
                    while (swap == i)
                        swap = rnd.Next(0, Cities.Length);
                    temp = perm[i];
                    perm[i] = perm[swap];
                    perm[swap] = temp;
                }
                Route.Clear();
                for (i = 0; i < Cities.Length; i++)                            // Now build the route using the random permutation 
                {
                    Route.Add(Cities[perm[i]]);
                }
                bssf = new TSPSolution(Route);
                count++;
            } while (costOfBssf() == double.PositiveInfinity);                // until a valid route is found
            timer.Stop();

            results[COST] = costOfBssf().ToString();                          // load results array
            results[TIME] = timer.Elapsed.ToString();
            results[COUNT] = count.ToString();

            return results;
        }

        /// <summary>
        /// performs a Branch and Bound search of the state space of partial tours
        /// stops when time limit expires and uses BSSF as solution
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] bBSolveProblem()
        {
            //This is used to give proportional advantages to the level of a state so that it would come off the queue before ones with lower levels.
            int priorityFactor = 10 * Cities.Length;
            int statesCreated = 0;
            int statesPruned = 0;
            int maxStoredStates = 0;
            int numBssfUpdates = 0;
            string[] results = new string[3];
            Route = new ArrayList();
            int length = Cities.Length;
            //This is the total cost matrix of all the cities to all of its connecting cities.
            double[,] costs = new double [length,length];
            Stopwatch stopWatch = new Stopwatch();
            //Binary min heap
            BinaryHeap queue = new BinaryHeap();
            //This is the upperbound of the BSSF
            double upperBound = double.PositiveInfinity;

            stopWatch.Start();

            /*
             * Iterate through the cities and fill up the total cost matrix.
             * Time comp: O(n^2)
             * Space comp: O(n^2)
             */
            for (int row = 0; row < length; row++)
            {
                for (int col = 0; col < length; col++)
                {
                    if (row == col) costs[row, col] = double.PositiveInfinity;
                    else costs[row, col] = Cities[row].costToGetTo(Cities[col]);
                }
            }

            /* Do greedy algorithm starting from every city and save the best to get BSSF.
             * BestRoute will save the route of the best greedy solution and upperBound will save the best upperbound.
             * Time comp: O(n^3) because we do the greedy algorithm (O(n^2)) for every city
             * Space comp: O(n^2) because we only have to store the data for each iteration of the greedy algorithm
             */
            ArrayList bestRoute = new ArrayList();
            for(int city = 0; city < length; city++)
            {
                ArrayList currRoute = new ArrayList();
                currRoute.Add(city);
                int startPosition = city;
                double currGreedyUpperbound = getGreedyUpperbound(city, city, costs, currRoute, length, 0);
                if(currGreedyUpperbound < upperBound)
                {
                    upperBound = currGreedyUpperbound;
                    bestRoute.Clear();
                    bestRoute = currRoute;
                }
            }

            //count = how many solutions we found.
            int count = 1;
            
            //First we make the starting state and add it to the priority queue.
            int[] state1_tour = new int[1];
            state1_tour[0] = 0;
            State state0 = new State(0, 1, costs, 0, state1_tour, length);
            statesCreated++;
            queue.AddToAllNodes(state0);
            
            /* Add the first state to the priority queue.
             * Time comp: O(nlogn)
             */
            queue.Add(state0, state0.lowerBound - priorityFactor * state0.level * state0.level);

            /* Go through the best greedy route and update the bssf object so that we have an initial BSSF to compare the starting state's children to.
             * Time comp: O(n) because we go through it once for every city
             * Space comp: O(n) can't get more than the number of cities
             */ 
            foreach (int i in bestRoute)
            {
                Route.Add(Cities[i]);
            }
            bssf = new TSPSolution(Route);

            //If the greedy solution's upperbound is equal to the first state's lowerbound, we know that the calculated greedy solution is the optimal solution.
            if (upperBound == state0.lowerBound)
            {
                Console.WriteLine("States Created: " + statesCreated);
                Console.WriteLine("States Pruned: " + statesPruned);
                Console.WriteLine("Max # of Stored States at a given time: " + maxStoredStates);
                Console.WriteLine("Times BSSF Updated: " + numBssfUpdates);
                Console.WriteLine("Priority Factor: " + priorityFactor);
                results[COST] = costOfBssf().ToString();
                results[TIME] = stopWatch.Elapsed.ToString();
                results[COUNT] = count.ToString();
                return results;
            }

            //while loop that stops once the time is up or once we have been through the whole state tree
            // Time complexity of the branch is O(n^2 * (n-1)!)

            while(stopWatch.ElapsedMilliseconds < time_limit && queue.Count > 0)
            {
                //store the max number of stored states at a given time
                if (queue.Count > maxStoredStates) maxStoredStates = queue.Count;

                /* Grabs the smallest priority state from the queue.
                 * Time comp: O(nlogn)
                 */
                State state = queue.ExtractMin();

                //If the BSSF solution's upperbound is equal to the current state's lowerbound, we know that the calculated BSSF solution is the optimal solution.
                if (state.lowerBound == upperBound)
                {
                    Console.WriteLine("States Created: " + statesCreated);
                    Console.WriteLine("States Pruned: " + statesPruned);
                    Console.WriteLine("Max # of Stored States at a given time: " + maxStoredStates);
                    Console.WriteLine("Times BSSF Updated: " + numBssfUpdates);
                    Console.WriteLine("Priority Factor: " + priorityFactor);
                    results[COST] = costOfBssf().ToString();
                    results[TIME] = stopWatch.Elapsed.ToString();
                    results[COUNT] = count.ToString();
                    return results;
                }

                //If the BSSF solution's upperbound is lower than the current state's lowerbound, we know that we can prune it and not expand its children.
                else if (state.lowerBound > upperBound) statesPruned++;
                else //Let's get the state's children!
                {
                    //If we're at a leaf and its bound is less than the current best so far bound, then we make this solution the new BSSF.
                    if (state.tour.Length == length + 1)
                    {
                        count++;
                        if (state.lowerBound < upperBound)
                        {
                            /* Go through the best route so far and update the bssf object.
                             * Time comp: O(n) because we go through it once for every city
                             * Space comp: O(n) can't get more than the number of cities
                             */ 
                            Route = new ArrayList();
                            for (int i = 0; i < state.tour.Length; i++)
                            {
                                Route.Add(Cities[state.tour[i]]);
                            }
                            bssf = new TSPSolution(Route);
                            numBssfUpdates++;
                            upperBound = state.lowerBound;
                        }
                    }
                    else //Get children and add them to the queue if they are not prunable.
                    {
                        /* Go through all of the current state's children while also getting the tour for each child.
                         * Time comp: O(n^2) because we calculate the tour for each children for the state
                         * Space comp: O(n^2) because we store all the children and their tour
                         */
                        for (int col = 0; col < length; col++)
                        {
                            double currCost = costs[state.pointsIndex, col];
                            if (!currCost.Equals(double.PositiveInfinity))
                            {
                                //Get the tour for the child.
                                int[] currTour = new int[state.level + 1];
                                for (int i = 0; i < state.tour.Length; i++)
                                {
                                    currTour[i] = state.tour[i];
                                }
                                currTour[state.level] = col;

                                //Build a state for the child.
                                State childState = new State(col, state.level + 1, state.cost, state.lowerBound, currTour, length);
                                statesCreated++;
                                queue.AddToAllNodes(childState);

                                //If it has potential to become the BSSF, add it to the queue so that we can expand it.
                                if (childState.lowerBound < upperBound)
                                {
                                    /* Add the first state to the priority queue.
                                     * Time comp: O(nlogn)
                                     */
                                    queue.Add(childState, childState.lowerBound - priorityFactor * childState.level * childState.level);
                                }

                                //If the BSSF solution's upperbound is lower than the current child's lowerbound, we know that we can prune it and not expand it.
                                else if (childState.lowerBound > upperBound) statesPruned++;
                            }
                        }
                    }
                }
            }
            
            stopWatch.Stop();

            //Either we ran out of time or we definitely have the optimal solution. Return the BSSF results.
            Console.WriteLine("States Created: " + statesCreated);
            Console.WriteLine("States Pruned: " + statesPruned);
            Console.WriteLine("Max # of Stored States at a given time: " + maxStoredStates);
            Console.WriteLine("Times BSSF Updated: " + numBssfUpdates);
            Console.WriteLine("Priority Factor: " + priorityFactor);
            results[COST] = costOfBssf().ToString();    // load results into array here, replacing these dummy values
            results[TIME] = stopWatch.Elapsed.ToString();
            results[COUNT] = count.ToString();
            return results;
        }

        /* The recursive greedy algorithm that travels through the graph and takes the greedy short path at each step.
         * Time comp: O(n^2) because we are checking each city's children for the greedy path after we take that path, we do it again for that city.
         * Space comp: O(n^2) because we store the cost matrix everytime.
         */
        private double getGreedyUpperbound(int startPosition, int currRow, double[,] costs, ArrayList currRoute, int numCities, double routeCost)
        {
            double minTo = double.PositiveInfinity;
            int to = -1;
            for (int i = 0; i < numCities; i++)
            {
                if (costs[currRow, i] < minTo && !currRoute.Contains(i))
                {
                    minTo = costs[currRow, i];
                    to = i;
                }
            }

            if (to == -1)
            {
                if (costs[currRow, startPosition] < double.PositiveInfinity && currRoute.Count == numCities)
                {
                    routeCost += costs[currRow, startPosition];
                    return routeCost;
                }
                return double.PositiveInfinity;
            }

            currRoute.Add(to);
            routeCost += minTo;
            currRow = to;
            return getGreedyUpperbound(startPosition, currRow, costs, currRoute, numCities, routeCost);
        }

        /////////////////////////////////////////////////////////////////////////////////////////////
        // These additional solver methods will be implemented as part of the group project.
        ////////////////////////////////////////////////////////////////////////////////////////////

        /// <summary>
        /// finds the greedy tour starting from each city and keeps the best (valid) one
        /// </summary>
        /// <returns>results array for GUI that contains three ints: cost of solution, time spent to find solution, number of solutions found during search (not counting initial BSSF estimate)</returns>
        public string[] greedySolveProblem()
        {
            string[] results = new string[3];
            Route = new ArrayList();
            int length = Cities.Length;
            //This is the total cost matrix of all the cities to all of its connecting cities.
            double[,] costs = new double[length, length];
            Stopwatch stopWatch = new Stopwatch();
            //Binary min heap
            BinaryHeap queue = new BinaryHeap();
            //This is the upperbound of the BSSF
            double upperBound = double.PositiveInfinity;

            stopWatch.Start();

            /*
             * Iterate through the cities and fill up the total cost matrix.
             * Time comp: O(n^2)
             * Space comp: O(n^2)
             */
            for (int row = 0; row < length; row++)
            {
                for (int col = 0; col < length; col++)
                {
                    if (row == col) costs[row, col] = double.PositiveInfinity;
                    else costs[row, col] = Cities[row].costToGetTo(Cities[col]);
                }
            }

            /* Do greedy algorithm starting from every city and save the best to get BSSF.
             * BestRoute will save the route of the best greedy solution and upperBound will save the best upperbound.
             * Time comp: O(n^3) because we do the greedy algorithm (O(n^2)) for every city
             * Space comp: O(n^2) because we only have to store the data for each iteration of the greedy algorithm
             */
            ArrayList bestRoute = new ArrayList();
            for (int city = 0; city < length; city++)
            {
                ArrayList currRoute = new ArrayList();
                currRoute.Add(city);
                int startPosition = city;
                double currGreedyUpperbound = getGreedyUpperbound(city, city, costs, currRoute, length, 0);
                if (currGreedyUpperbound < upperBound)
                {
                    upperBound = currGreedyUpperbound;
                    bestRoute.Clear();
                    bestRoute = currRoute;
                }
            }

            /* Go through the best greedy route and update the bssf object so that we have an initial BSSF to compare the starting state's children to.
             * Time comp: O(n) because we go through it once for every city
             * Space comp: O(n) can't get more than the number of cities
             */ 
            foreach (int i in bestRoute)
            {
                Route.Add(Cities[i]);
            }

            stopWatch.Stop();
            bssf = new TSPSolution(Route);
            results[COST] = costOfBssf().ToString();    // load results into array here, replacing these dummy values
            results[TIME] = stopWatch.Elapsed.ToString();
            results[COUNT] = "1"; //*** should this be 1 or the number of greedy solutions that we find???

            return results;
        }

        public string[] fancySolveProblem()
        {
            greedySolveProblem();
            return new GeneticAlgorithm().Solve(this);
        }
        #endregion
    }
    /*
     * State class for use in the binary heap.
     * Parameters: index of the city, which level in the branch and bound tree we are at, 
     *              the parent's cost matrix, the parent's lowerbound, the tour of this current state, and total number of cities total.
     */
    public class State
    {
        public int pointsIndex, QueuePosition, level;
        public double distance, lowerBound;
        public double[,] cost;
        public int[] tour;

        /* Building a state. 
         * Time comp: O(n^2) because we reduce the parent's cost matrix by iterating over each row and col
         * Space comp: O(n^2) because we store the reduced cost matrix within the state 
         */
        public State(int ptIndex, int level, double[,] parentCosts, double parentLowerbound, int[] tour, int numCities)
        {
            this.level = level;
            this.tour = tour;
            this.pointsIndex = ptIndex;
            this.QueuePosition = -1;
            this.distance = 0;

            //If it is the first state made, the lowerbound is 0.
            if (level == 1) this.lowerBound = 0;
            else // Otherwise it is equal to the parent's lowerbound + the cost of the edge to get to this city + future changes to the reduced matrix.
            {
                this.lowerBound = parentLowerbound + parentCosts[tour[tour.Length - 2], ptIndex];
            }

            //Copy the parent's cost matrix O(n^2)
            this.cost = new double[numCities, numCities];
            for(int row = 0; row < numCities; row++)
                for(int col = 0; col < numCities; col++)
                    this.cost[row, col] = parentCosts[row, col];

            /* We need to put infinities on the row of the city we came from and the column of the city that we went to. Also on the reflexive cell. (cost[col][col])
             * Time comp: O(n)
             */
            if(level > 1)
            {
                int cameFrom = tour[tour.Length - 2];
                for(int col = 0; col < numCities; col++)
                    this.cost[cameFrom, col] = double.PositiveInfinity;
                for(int row = 0; row < numCities; row++)
                    this.cost[row, ptIndex] = double.PositiveInfinity;
                this.cost[ptIndex, cameFrom] = double.PositiveInfinity;
            }

            /* We go through the cost matrix row by row and find the minimum for each row and save the differences (cell - row minimum) in each of the row's cells.
             * Time comp: O(n^2) because we check and change every value in the 2D matrix
             */
            for(int row = 0; row < numCities; row++)
            {
                double rowMin = double.PositiveInfinity;
                for(int col = 0; col < numCities; col++)
                {
                    double currCost = this.cost[row, col];
                    if (!currCost.Equals(double.PositiveInfinity) && currCost < rowMin)
                    {
                        rowMin = currCost;
                    }
                }
                for (int col = 0; col < numCities; col++)
                {
                    if (!this.cost[row, col].Equals(double.PositiveInfinity)) this.cost[row, col] = cost[row, col] - rowMin;
                }
                if (rowMin < double.PositiveInfinity)
                {
                    this.lowerBound += rowMin;
                }
            }

            /* We go through the cost matrix column by column and find the minimum for each column and save the differences (cell - row minimum) in each of the column's cells.
             * We do this if there is no zero in a column left from the previous step. This reduces the cost matrix more.
             * Time comp: O(n^2) because we check and possibly change every value in the 2D matrix
             */
            for (int col = 0; col < numCities; col++)
            {
                bool isZero = false;
                double colMin = double.PositiveInfinity;
                for (int row = 0; row < numCities; row++)
                {
                    double currCost = cost[row, col];
                    if (currCost == 0)
                    {
                        isZero = true;
                        break;
                    }
                    else if (currCost < colMin) colMin = currCost;
                }
                if(!isZero)
                {
                    for (int row = 0; row < numCities; row++)
                    {
                        if (!this.cost[row, col].Equals(double.PositiveInfinity)) this.cost[row, col] = this.cost[row, col] - colMin;
                    }
                    if (colMin < double.PositiveInfinity)
                    {
                        this.lowerBound += colMin;
                    }
                }
            }
        }
    }
    /// <summary>
    /// A min-type priority queue of Nodes
    /// </summary>
    public class BinaryHeap
    {
        #region Instance variables

        public State[] allNodes;
        readonly State[] data;
        readonly double[] distances;
        public int count;
        #endregion

        /// <summary>
        /// Creates a new, empty priority queue with the specified capacity.
        /// </summary>
        /// <param name="capacity">The maximum number of nodes that will be stored in the queue.</param>
        public BinaryHeap()
        {
            //Have this be our capacity for the heap. If we do less than 50 cities, this should be enough space.
            int capacity = 10000000;
            allNodes = new State[capacity];
            data = new State[capacity];
            distances = new double[capacity];
            count = 0;
        }

        /// <summary>
        /// Adds an item to the queue.  Is position is determined by its priority relative to the other items in the queue.
        /// aka HeapInsert
        /// </summary>
        /// <param name="item">Item to add</param>
        /// <param name="priority">Priority value to attach to this item.  Note: this is a min heap, so lower priority values come out first.</param>
        public void Add(State item, double priority)
        {
            if (count == data.Length)
                throw new Exception("Heap capacity exceeded");

            // Add the item to the heap in the end position of the array (i.e. as a leaf of the tree)
            int position = count++;
            data[position] = item;
            item.QueuePosition = position;
            distances[position] = priority;
            // Move it upward into position, if necessary
            MoveUp(position);
        }

        /// <summary>
        /// Add a state to the heap's copy of the total made states so that we have a list of every state made.
        /// </summary>
        /// <param name="item">Item to add</param>
        public void AddToAllNodes(State item)
        {
            allNodes[item.pointsIndex] = item;
        }

        /// <summary>
        /// Extracts the item in the queue with the minimal priority value.
        /// </summary>
        /// <returns></returns>
        public State ExtractMin()
        {
            allNodes[data[0].pointsIndex].distance = distances[data[0].QueuePosition];
            State minNode = data[0];
            Swap(0, count - 1);
            count--;
            MoveDown(0);
            return minNode;
        }

        /// <summary>
        /// Reduces the priority of a node already in the queue.
        /// aka DecreaseKey 
        /// </summary>
        public void DecreasePriority(int index, double priority)
        {
            int position = allNodes[index].QueuePosition;
            while ((position > 0) && (distances[Parent(position)] > priority))
            {
                int original_parent_pos = Parent(position);
                Swap(original_parent_pos, position);
                position = original_parent_pos;
            }
            distances[position] = priority;
        }

        /// <summary>
        /// Moves the node at the specified position upward, it it violates the Heap Property.
        /// This is the while loop from the HeapInsert procedure in the slides.
        /// </summary>
        /// <param name="position"></param>
        void MoveUp(int position)
        {
            while ((position > 0) && (distances[Parent(position)] > distances[position]))
            {
                int original_parent_pos = Parent(position);
                Swap(position, original_parent_pos);
                position = original_parent_pos;
            }
        }

        /// <summary>
        /// Moves the node at the specified position down, if it violates the Heap Property
        /// aka Heapify
        /// </summary>
        /// <param name="position"></param>
        void MoveDown(int position)
        {
            int lchild = LeftChild(position);
            int rchild = RightChild(position);
            int largest = 0;
            if ((lchild < count) && (distances[lchild] < distances[position]))
            {
                largest = lchild;
            }
            else
            {
                largest = position;
            }
            if ((rchild < count) && (distances[rchild] < distances[largest]))
            {
                largest = rchild;
            }
            if (largest != position)
            {
                Swap(position, largest);
                MoveDown(largest);
            }
        }

        /// <summary>
        /// Number of items waiting in queue
        /// </summary>
        public int Count
        {
            get
            {
                return count;
            }
        }

        #region Utilities
        /// <summary>
        /// Swaps the nodes at the respective positions in the heap
        /// Updates the nodes' QueuePosition properties accordingly.
        /// </summary>
        void Swap(int position1, int position2)
        {
            State temp = data[position1];
            data[position1] = data[position2];
            data[position2] = temp;
            data[position1].QueuePosition = position1;
            data[position2].QueuePosition = position2;

            double temp2 = distances[position1];
            distances[position1] = distances[position2];
            distances[position2] = temp2;
        }

        /// <summary>
        /// Gives the position of a node's parent, the node's position in the queue.
        /// </summary>
        static int Parent(int position)
        {
            return (position - 1) / 2;
        }

        /// <summary>
        /// Returns the position of a node's left child, given the node's position.
        /// </summary>
        static int LeftChild(int position)
        {
            return 2 * position + 1;
        }

        /// <summary>
        /// Returns the position of a node's right child, given the node's position.
        /// </summary>
        static int RightChild(int position)
        {
            return 2 * position + 2;
        }

        /// <summary>
        /// Checks all entries in the heap to see if they satisfy the heap property.
        /// </summary>
        public void TestHeapValidity()
        {
            for (int i = 1; i < count; i++)
                if (distances[Parent(i)] > distances[i])
                    throw new Exception("Heap violates the Heap Property at position " + i);
        }
        #endregion
    }
}
