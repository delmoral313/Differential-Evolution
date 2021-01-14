using System;
class DifferentialEvolution
{
    //Declare all variables and data types needed for the bioinspired algorithm
    ushort popsize,dim;
    float Cr,F,lb,ub;
    struct individual{public float[] x; public float fitness;}; //Custom data type
    individual[] population,next_population;
    Random rnd;
    public delegate float ObjectiveFunction(float[] x);
    public DifferentialEvolution(ushort dim,ushort popsize,ushort seed,float Cr,float F,float lb,float ub)
    {
        //Constructor 
        this.popsize = popsize;
        this.dim = dim;
        this.rnd = new Random(seed);
        this.Cr = Cr;
        this.F = F;
        this.lb = lb; //Lower bound
        this.ub = ub; //Upper bound
    }
    public void InitializePopulation(ObjectiveFunction fobj)
    {
        population = new individual[this.popsize];
        next_population = new individual[this.popsize];
        for (ushort i = 0; i < this.popsize; i++)
        {
            population[i].x = new float[this.dim];
            next_population[i].x = new float[this.dim];
            for (ushort j = 0; j < this.dim; j++)
            {
                //Generate a number in range [lb,ub]
                population[i].x[j] = (float)(this.lb + rnd.NextDouble()*(this.ub-this.lb));
            }
            /*
            -1 since we are minimizing the fobj but maximizing fitness. If we were maximizing 
            fobj there'd be no need for the -1 factor.  
            */
            population[i].fitness = (-1)*fobj(population[i].x); //greater is better 
            //Console.WriteLine("{0}:{1}",i,population[i].fitness);
        }
    }
    public void EvolvePopulation(ObjectiveFunction fobj,uint maxGen)
    {
        //maxFEs is maximum number of function evaluations 
        individual trial; 
        uint G = 0;
        trial.x = new float[this.dim]; 
        while(G < maxGen)
        {
            for (ushort i = 0; i < this.popsize; i++) //for every individual 
            {   
                this.Mutate(trial.x,i);
                this.CrossOver(trial.x,i);
                trial.fitness = (-1)*fobj(trial.x); 
                if(trial.fitness >= population[i].fitness) 
                {
                    next_population[i].fitness = trial.fitness;
                    trial.x.CopyTo(next_population[i].x,0);
                }
                else 
                {
                    next_population[i].fitness = population[i].fitness;
                    population[i].x.CopyTo(next_population[i].x,0);
                }
            }
            //Advance generation
            this.next_population.CopyTo(this.population,0);
            G += 1;
        }
    }
    public void GetFittest()
    {
        //Find fittest candidate solution in O(n)
        float max_fitness = this.population[0].fitness;
        ushort candidate_id = 0;
        float sum_fitness = this.population[0].fitness;
        for (ushort i = 1; i < this.popsize; i++)
        {
            if(this.population[i].fitness > max_fitness) 
            {
                max_fitness = this.population[i].fitness;
                candidate_id = i;
            }
            sum_fitness += this.population[i].fitness;
        }
        Console.WriteLine("Candidate Id = {0}\tf(xo) = {1}",candidate_id,-max_fitness);
        Console.WriteLine("Mean f(x)= {0}",-sum_fitness/this.popsize);
    }
    void Mutate(float[] v, ushort i)
    {
        //Force diversity
        ushort a = (ushort)this.rnd.Next(this.popsize); while(a == i) a = (ushort)this.rnd.Next(this.popsize);
        ushort b = (ushort)this.rnd.Next(this.popsize); while(b == a || b == i) b = (ushort)this.rnd.Next(this.popsize);
        ushort c = (ushort)this.rnd.Next(this.popsize); while (c == a || c == b || c == i) c = (ushort)this.rnd.Next(this.popsize);
        //Vector operation: u = pop[a] + F*(pop[b]-pop[c])
        for (ushort j = 0; j < this.dim; j++)
        {
            v[j] = this.population[a].x[j] + this.F*(this.population[b].x[j] - this.population[c].x[j]);
            //Take care of bounds 
            if(v[j] < this.lb)
                v[j] = (this.lb + this.population[i].x[j])/2;
            else if(v[j] > this.ub)
                v[j] = (this.ub + this.population[i].x[j])/2;
        }
    }
    void CrossOver(float[] u,ushort i)
    {   //Binomial recombination 
        ushort jrand = (ushort)this.rnd.Next(this.dim);
        for (ushort j = 0; j < this.dim; j++)
        {
            if(this.rnd.NextDouble()<=this.Cr || j==jrand) continue;
            else u[j] = this.population[i].x[j]; 
        }
    }
}

class Program
{
    static void Main(string[] args)
    {
        /*
        the only way to calling a non-static method within the same class 
        Program prog = new Program(); 
        prog.print_something();
        */
        ushort dim = 100;
        ushort popsize = 200;
        ushort seed = 7;
        float Cr = (float)0.9;
        float F = (float)0.4;
        float lb = -100;
        float ub = 100;
        uint maxGen = 4000;
        DifferentialEvolution DE = new DifferentialEvolution(dim,popsize,seed,Cr,F,lb,ub);
        DE.InitializePopulation(my_ObjectiveFunction);
        DE.EvolvePopulation(my_ObjectiveFunction,maxGen);
        DE.GetFittest();
    }
    static float my_ObjectiveFunction(float[] x)
    {
        //F12 from CEC'13 LSGO Benchmark 
        float sum = 0;
        float term1,term2;
        for (ushort j = 0; j < x.Length-1; j++)
        {
            term1 = ((x[j]*x[j])-x[j+1])*((x[j]*x[j])-x[j+1]);
            term2 = (x[j]-1)*(x[j]-1);
            sum += (100*term1 + term2);
        }
        return sum;
    }
}
