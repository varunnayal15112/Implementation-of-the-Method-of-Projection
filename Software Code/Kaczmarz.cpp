	/*
		Program to implement Method Of Projection for Algebraic Reconstruction Algorithms.
		Using Kaczmarz Algorithm/Method.
		Filename : Kaczmarz.cpp
		Date of last update : 09 Jan 2018
	*/


/*****HEADER FILES*****/
//#include<bits/stdc++.h>
#include<iostream.h>
#include <time.h>

/*****ASSUMPTIONS*****/
typedef long int ll;
#define no_iteration 100000                                // MAXIMUM LIMIT OF THE ITERATIONS 
#define MAX 1000										   // MAXIMUM LIMIT OF THE RAYS(EQUATIONS)
//using namespace std;

/*****VARIABLE DECLARATION*****/
ll i,j,k;												   // VARIABLES FOR CONTROLLING LOOPS
ll no_variable,no_equation,row_of_convergence=no_iteration;

	/*
		no_variable        -> number of unknowns(gray levels)
		no_equation        -> number of rays
		row_of_convergence -> row where convergence achieved,initially initialized with the maximum iterations
	*/

/*****INITIALIZE VARIABLES AS PER USER CHOICE*****/
void get_variable_from_user()
{
    cout<<"\nEnter number of Rays/Equations : ";
    cin>>no_equation;
    cout<<"\nEnter number of unknowns/gray_levels : ";
    cin>>no_variable;
}

/*****DATA STRUCTURE DECLARATION*****/
long double weight_factor_array[MAX][MAX];
long double ray_sum_array[MAX];
long double pixel_gray_value[no_iteration+1][MAX];
long double previous_pixel_gray_value[MAX];
long double unknown[MAX];

    /*
		weight_factor_array[][]     -> 2-D array to store the weighting-factor of the rays  
		ray_sum_array[]             -> 1-D array to store the ray_sums of the rays
		pixel_gray_value[][]        -> 2-D array to store the gray levels at each iteration
		previous_pixel_gray_value[] -> 1-D array to store the gray levels of previos projection
		unknown[]                   -> 1-D array to store the grey levels after each iteration
    */

/*****GETTING DATA FROM USER LIKE WEIGHTING FACTOR AND RAY-SUM OF THE RAYS*****/
void get_equation_from_user()
{
    cout<<"\nEnter the weighting-factor of the equations : ";
    for(i=1;i<=no_equation;i++)
    {
        cout<<"\nEquation number "<<i<<" -> ";
        for(j=0;j<no_variable;j++)
            cin>>weight_factor_array[i][j];
    }
    cout<<"\nEnter the Ray-Sum of the equations : ";
    for(i=1;i<=no_equation;i++)
    {
        cout<<"\nEquation number "<<i<<" -> ";
        cin>>ray_sum_array[i];
    }
}

/*****INITIALIZE THE INITIAL GUESS POINT WITH ALL ZEROES*****/
void initialize_initial_pixelGray()
{
    for(i=0;i<no_variable;i++)
    {
        pixel_gray_value[0][i]=0.0;
    }
}

/*****PRINT USER DEFINED DATA ON CONSOLE*****/
void print_details()
{
    cout<<"\nNumber of Rays/Equations       : "<<no_equation;
    cout<<"\nNumber of Unknowns/gray_levels : "<<no_variable;
    cout<<"\nMaximum Limit of iterations    : "<<no_iteration;
    cout<<"\n...................Information of Rays...................\n";
	for(i=1;i<=no_equation;i++)
    {
		cout<<"\nWeighting-factor of Ray/Equation "<<i<<" : ";
        for(j=0;j<no_variable;j++)
            cout<<weight_factor_array[i][j]<<" ";
    }
    for(i=1;i<=no_equation;i++)
		cout<<"\nRay-Sum of Ray/Equation "<<i<<" : "<<ray_sum_array[i];
}

/*****PROCEDURES TO CALCULATE THE VALUES OF THE KACZMARZ FORMULAE*****/
long double calculate_numerator(ll eq_no)
{
    long double sum=0.0;
    int x;
    for(x=0;x<no_variable;x++)
    {
        sum=sum+(previous_pixel_gray_value[x]*weight_factor_array[eq_no][x]);
    }
    sum=sum-ray_sum_array[eq_no];
    return sum;
}

long double calculate_denominator(ll eq_no)
{
    long double sum=0.0;
    int x;
    for(x=0;x<no_variable;x++)
    {
        sum=sum+weight_factor_array[eq_no][x]*weight_factor_array[eq_no][x];
    }
    return sum;
}

long double calculate_multiplier(ll eq_no,ll var_no)
{
    return (weight_factor_array[eq_no][var_no]);
}

/*****KACZMARZ ALGORITHM IMPLEMENTATION*****/
void implement_kaczmarz_method()
{
    for(i=0;i<no_variable;i++)
    {
        previous_pixel_gray_value[i]=pixel_gray_value[0][i];
    }
    for(i=1;i<=no_iteration;i++)
    {
        for(j=1;j<=no_equation;j++)
        {
            for(k=0;k<no_variable;k++)
            {
                unknown[k]=previous_pixel_gray_value[k]-((calculate_numerator(j)/calculate_denominator(j))*calculate_multiplier(j,k));
            }
            for(k=0;k<no_variable;k++)
            {
                previous_pixel_gray_value[k]=unknown[k];
            }
        }
        for(k=0;k<no_variable;k++)
        {
            pixel_gray_value[i][k]=previous_pixel_gray_value[k];
        }
        int flag=0;
        for(k=0;k<no_variable;k++)
        {
            if(pixel_gray_value[i][k]!=pixel_gray_value[i-1][k])
            {
                flag=1;
                break;
            }
        }
        if(!flag)
        {
			//convergence achieved before the maximum limit of iteration
            row_of_convergence=i;
            break;
        }
    }
    
    /*
    for(i=0;i<=no_iteration;i++)
    {
        for(k=0;k<no_variable;k++)
            cout<<pixel_gray_value[i][k]<<" ";
        cout<<"\n";
    }
    */   
}

/*****PRINT SOLUTION ON CONSOLE*****/
void put_solution_on_console()
{
    cout<<"\nSolution of given rays/equations is : \n";
    for(i=0;i<no_variable;i++)
    {
        cout<<pixel_gray_value[row_of_convergence][i]<<" ";
    }
    cout<<"\nIterations Required : "<<row_of_convergence;
}

/*****MAIN FUNCTION -> ENTRY POINT OF THE PROGRAM*****/
int main()
{
	
	clock_t t1, t2;

    get_variable_from_user();
    get_equation_from_user();
    
	//starting timer
	t1 = clock();
	
	initialize_initial_pixelGray();
    print_details();
    implement_kaczmarz_method();
    put_solution_on_console();

	//stop timer
	t2 = clock();
	cout<<"\n\nTime taken in computations : "<<((double)(t2 - t1)/CLOCKS_PER_SEC)<<" seconds\n\n";

    return 0;
}