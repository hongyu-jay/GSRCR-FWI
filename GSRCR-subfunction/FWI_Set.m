
function  [ Opts]   =    FWI_Set(   I, name, lambda, alpha, beta, mu,  c1, c2)

randn ('seed',0);

Opts.I         =   double(I);

Opts.Iter      =   35;

Opts.step      =   4;         

Opts.patch     =  7;   

Opts.Sim       =  40;  


if name == 1
      
    Opts.lambda      =  lambda;   
    
    Opts. alpha      =  alpha;    
    
    Opts. beta       =  beta; 
    
    Opts. mu         =  mu;     
    
    Opts. c1         =   c1;   

    Opts. c2         =   c2; 
    Opts.IterNums = 20;
Opts.patch = 7;
Opts.Region = 25;
Opts.Sim = 60;
Opts.step = 4;
Opts.paraell = 1;
    
elseif name == 2
    
    Opts.lambda      =  lambda;   
    
    Opts. alpha      =  alpha;    
    
    Opts. beta       =  beta; 
    
    Opts. mu         =  mu;     
    
    Opts. c1         =   c1;   

    Opts. c2         =   c2;    
    
    Opts.alpha = alpha;
Opts.IterNums = 20;
Opts.patch = 6;
Opts.Region = 35;
Opts.Sim = 40;
Opts.step = 1;
Opts.paraell = 1;   
     
end


end

