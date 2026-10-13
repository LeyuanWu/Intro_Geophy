##################################################################
# # ! Leyuan Wu
# Functions to do statistics
##################################################################

import numpy as np;

def DispVars(Vars, Names, Units):
    """Display statistics of a list of np.array variables

    Args:
        Vars (numpy array list) : a list of np.array variables
        Names (str list)        : Physical quantity names
        Units (str list)        : Units
    """
    for var,name,unit in zip(Vars,Names,Units):
        print('\n*********************',flush=True);
        print('%s'%(name),flush=True);
        print('*********************',flush=True);
        if var.dtype != np.dtype('datetime64[s]'):
            tp_var = type(var);
            shp_var = np.shape(var);
            min_var = np.min(var); 
            max_var = np.max(var); 
            mean_var = np.mean(var); 
            STD = np.sqrt(np.mean((var-mean_var)**2));
            print(f'<Data type>: {tp_var}', flush=True);
            print(f'<Shape>: {shp_var}', flush=True);
            print('<Min>:  %12.4f [%s]'%(min_var,unit),flush=True);
            print('<Max>:  %12.4f [%s]'%(max_var,unit),flush=True);
            print('<Mean>: %12.4f [%s]'%(mean_var,unit),flush=True);
            print('<STD>:  %12.4f [%s]'%(STD,unit),flush=True);
        else:
            tp_var = type(var);
            shp_var = np.shape(var);
            min_var = np.min(var); 
            max_var = np.max(var); 
            print(f'<Data type>: {tp_var}', flush=True);
            print(f'<Shape>: {shp_var}', flush=True);
            print(f'<Min>:  {min_var}', flush=True);
            print(f'<Max>:  {max_var}', flush=True);
        print('*********************\n',flush=True);
    
def DispSTATS(Est,Ref,Name,unit):
    """Display statistics of the difference between two data

    Args:
        Est (numpy array): Estimated value
        Ref (numpy array): Reference value
        Name (str): Physical quantity
        unit (str): Unit
    """
    dif = Est-Ref;
    minErr = np.min(dif); 
    maxErr = np.max(dif); 
    meanErr = np.mean(dif); 
    STD = np.sqrt(np.mean((dif-meanErr)**2));
    MSE   = np.mean(dif**2);
    RMSE  = np.sqrt(MSE);
    print('\n*********************',flush=True);
    print('%s'%(Name),flush=True);
    print('*********************',flush=True);
    print('<Min>:  %12.4f [%s]'%(minErr,unit),flush=True);
    print('<Max>:  %12.4f [%s]'%(maxErr,unit),flush=True);
    print('<Mean>: %12.4f [%s]'%(meanErr,unit),flush=True);
    print('<STD>:  %12.4f [%s]'%(STD,unit),flush=True);
    print('<MSE>:  %12.4f [%s]^2'%(MSE,unit),flush=True);
    print('<RMSE>: %12.4f [%s]'%(RMSE,unit),flush=True);
    L2NormRef = np.sqrt(np.mean((Ref)**2));
    if L2NormRef != 0:
        E2 = RMSE / L2NormRef;
        print('<E2>: %12.4e = %8.4f%%'%(E2,100*E2),flush=True);
    print('*********************\n',flush=True);
    