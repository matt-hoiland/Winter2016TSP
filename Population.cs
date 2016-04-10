using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TSP
{
    internal class Population
    {
        private List<Sequence> pop;
        private bool needSort = true;
        private Random random = new Random();

        public int Size
        {
            get
            {
                return pop.Count;
            }
        }

        public Population()
        {
            pop = new List<Sequence>();
        }

        public void Add(Sequence s)
        {
            pop.Add(s);
            needSort = true;
        }

        public Sequence this[int i]
        {
            get
            {
                return pop[i];
            }
        }

        public Sequence Best
        {
            get
            {
                if (needSort)
                {
                    pop.Sort();
                    needSort = false;
                }
                return pop[0];
            }
        }

        public Sequence GetCandidate()
        {
            return pop[random.Next(pop.Count)];
        }

        public void ForEach(Action<Sequence> doit)
        {
            foreach (var entry in pop)
            {
                doit(entry);
                needSort = true;
            }
        }

        public void TrimDuplicates()
        {
          
        }

        public void TrimWorst()
        {
            
        }
    }
}
