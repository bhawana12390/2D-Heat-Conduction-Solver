#include <iostream>
#include <fstream>
#include<bits/stdc++.h>
#define ll long long
using namespace std;

void apply_boundary_conditions(vector<vector<float>>& temp,float temp_x_0,float temp_x_length,float temp_y_0,float temp_y_height){
    for(int i=1;i<temp.size()-1;i++){
        temp[i][0]=temp_x_0;
        temp[i][temp[i].size()-1]=temp_x_length;
    }
    for(int i=1;i<temp[0].size()-1;i++){
        temp[0][i]=temp_y_0;
        temp[temp.size()-1][i]=temp_y_height;
    }
    for(int i=1;i<temp.size()-1;i++){
        for(int j=1;j<temp[i].size()-1;j++){
            temp[i][j]=0;
        }
    }
}
double calculate_residual(vector<vector<float>> temp,double beta){
    double res=0;
    for(int i=1;i<temp.size()-1;i++){
        for(int j=1;j<temp[0].size()-1;j++){
            double temp_temp=temp[i][j+1]+temp[i][j-1]+beta*temp[i+1][j]+beta*temp[i-1][j]-((float)(2*(1+beta))*temp[i][j]);
            res+=(temp_temp*temp_temp);
        }
    }
    return sqrt(res);
}
void solve_using_point_GS(vector<vector<float>>& temp,float over_relaxation,double beta,int max_iterations,double delta_x,double delta_y,vector<float>& residual_vector,vector<int>& comment,int j,float wanted_over_relaxation,int sym){
    // cout<<j<<" "<<comment[j]<<endl;
    int iterations=0;
    double residual = calculate_residual(temp,beta);
    // cout<<residual<<endl;
    if(over_relaxation<=(wanted_over_relaxation+0.01) && over_relaxation>=(wanted_over_relaxation-0.01)){
        residual_vector.push_back(residual);
    }
    while(iterations<max_iterations && residual>=min((delta_x*delta_x),(delta_y*delta_y))){
        iterations++;
        for(int i=1;i<temp.size()-1;i++){
            for(int j=1;j<temp[0].size()-1;j++){
                float temp_temp=temp[i][j+1]+temp[i][j-1]+beta*temp[i+1][j]+beta*temp[i-1][j]-((float)(2*(1+beta))*temp[i][j]);
                temp_temp*=(over_relaxation/(2*(1+beta)));
                temp[i][j]+=temp_temp;
            }
        }
        double prev_residual=residual;
        residual=calculate_residual(temp,beta);
        // cout<<residual<<" "<<iterations<<endl;
        if(over_relaxation<=(wanted_over_relaxation+0.01) && over_relaxation>=(wanted_over_relaxation-0.01)){
            residual_vector.push_back(residual);
        }
        if(prev_residual==residual || (prev_residual*1.1)<residual ){
            // cout<<prev_residual<<" "<<residual<<" "<<over_relaxation<<endl;
            comment[j]=1;
            break;
        }
        // for symmetric
        if(sym){
            iterations++;
            for(int i=temp.size()-2;i>=1;i--){
                for(int j=temp[0].size()-2;j>=1;j--){
                    float temp_temp=temp[i][j+1]+temp[i][j-1]+beta*temp[i+1][j]+beta*temp[i-1][j]-((float)(2*(1+beta))*temp[i][j]);
                    temp_temp*=(over_relaxation/(2*(1+beta)));
                    temp[i][j]+=temp_temp;
                }
            }
            prev_residual=residual;
            residual=calculate_residual(temp,beta);
            // cout<<residual<<" "<<iterations<<endl;
            if(over_relaxation<=(wanted_over_relaxation+0.01) && over_relaxation>=(wanted_over_relaxation-0.01)){
                residual_vector.push_back(residual);
            }
            if(prev_residual==residual || (prev_residual*1.1)<residual ){
                // cout<<prev_residual<<" "<<residual<<" "<<over_relaxation<<endl;
                comment[j]=1;
                break;
            }
        }
    }
    cout<<over_relaxation<<" "<<iterations<<" "<<comment[j]<<endl;
}
vector<vector<float>>set_coeff(vector<vector<float>> coefficients,vector<vector<float>>& temp,double beta){
    for(int i=1;i<coefficients.size()-1;i++){
        coefficients[i][0]=1;
        coefficients[i][1]=-2*(1+beta);
        coefficients[i][2]=1;
        // cout<<i<<" "<<coefficients[i][0]<<" "<<coefficients[i][1]<<" "<<coefficients[i][2]<<endl;
    }
    coefficients[0][0]=0;
    coefficients[0][1]=1;
    coefficients[0][2]=0;
    coefficients[coefficients.size()-1][0]=0;
    coefficients[coefficients.size()-1][1]=1;
    coefficients[coefficients.size()-1][2]=0;
    coefficients[1][0]=0;
    coefficients[coefficients.size()-2][2]=0;
    return coefficients;
}
vector<vector<float>>set_coeff_y(vector<vector<float>> coefficients,vector<vector<float>>& temp,double beta){
    for(int i=1;i<coefficients.size()-1;i++){
        coefficients[i][0]=beta;
        coefficients[i][1]=-2*(1+beta);
        coefficients[i][2]=beta;
        // cout<<i<<" "<<coefficients[i][0]<<" "<<coefficients[i][1]<<" "<<coefficients[i][2]<<endl;
    }
    coefficients[0][0]=0;
    coefficients[0][1]=1;
    coefficients[0][2]=0;
    coefficients[coefficients.size()-1][0]=0;
    coefficients[coefficients.size()-1][1]=1;
    coefficients[coefficients.size()-1][2]=0;
    coefficients[1][0]=0;
    coefficients[coefficients.size()-2][2]=0;
    return coefficients;
}
vector<float>set_rhs(vector<float>rhs,vector<vector<float>>& temp,double beta,int i){
    rhs[0]=temp[i][0];
    // cout<<"rhs"<<endl;
    rhs[rhs.size()-1]=temp[i][temp[i].size()-1];
    for(int j=1;j<temp[i].size()-1;j++){
        rhs[j]=-beta*(temp[i-1][j]+temp[i+1][j]);
        // cout<<beta<<" "<<temp[i-1][j]<<" "<<temp[i+1][j]<<" "<<rhs[j]<<endl;
    }
    rhs[1]-=temp[i][0];
    rhs[rhs.size()-2]-=temp[i][temp[i].size()-1];
    return rhs;
}
vector<float>set_rhs_y(vector<float>rhs,vector<vector<float>>& temp,double beta,int i){
    rhs[0]=temp[0][i];
    // cout<<"rhs"<<endl;
    rhs[rhs.size()-1]=temp[temp.size()-1][i];
    for(int j=1;j<rhs.size()-1;j++){
        rhs[j]=-(temp[j][i-1]+temp[j][i+1]);
        // cout<<beta<<" "<<temp[i-1][j]<<" "<<temp[i+1][j]<<" "<<rhs[j]<<endl;
    }
    rhs[1]-=(beta*temp[0][i]);
    rhs[rhs.size()-2]-=(beta*temp[temp.size()-1][i]);
    return rhs;
}
void apply_thomas_algorithm(vector<vector<float>>& coefficients, vector<float>&rhs){
    for(int i=2;i<coefficients.size()-1;i++){
        float multiplier= coefficients[i][0]/coefficients[i-1][1];
        coefficients[i][0]=0;
        coefficients[i][1]-=(coefficients[i-1][2]*multiplier);
        rhs[i]-=(rhs[i-1]*multiplier);
    }
    for(int i=coefficients.size()-3;i>=1;i--){
        float multiplier= coefficients[i][2]/coefficients[i+1][1];
        coefficients[i][2]=0;
        rhs[i]-=(rhs[i+1]*multiplier);
    }
}

void solve_using_line_GS(vector<vector<float>>& temp,double beta,int max_iterations,double delta_x,double delta_y,vector<float>& residual_vector_line,int sym,float wanted_over_relaxation,float over_relaxation,bool & line_gs){
    int iterations=0;
    double residual = calculate_residual(temp,beta); 
    double prev_residual=residual;
    int comment=0;
    int counter=0;
    if(abs(wanted_over_relaxation-over_relaxation)<=0.01){
        residual_vector_line.push_back(residual);
        // cout<<residual<<endl;
    }
    // cout<<iterations<<" "<<residual<<endl;
    vector<vector<float>>coefficients(temp[0].size(),vector<float>(3,0));
    vector<float>rhs(temp[0].size(),0);
    while(iterations<max_iterations && residual>=min((delta_x*delta_x),(delta_y*delta_y))){
        // cout<<"hehe"<<endl;
        iterations++;
        for(int i=1;i<temp.size()-1;i++){
            
            coefficients=set_coeff(coefficients,temp,beta);
            rhs=set_rhs(rhs,temp,beta,i);

            apply_thomas_algorithm(coefficients,rhs);
            for(int j=0;j<temp[i].size();j++){
                temp[i][j]+=(over_relaxation*(rhs[j]/coefficients[j][1]-temp[i][j]));
            }
            // cout<<endl;
        }
        residual=calculate_residual(temp,beta);
        if(abs(wanted_over_relaxation-over_relaxation)<=0.01){
            residual_vector_line.push_back(residual);
            // cout<<residual<<endl;
            if(prev_residual==residual || (prev_residual*1.1)<residual ){
                line_gs=false;
                comment=1;
                // cout<<"here";
                break;
            }
        }
        if(prev_residual==residual || (prev_residual*1.1)<residual){
            // cout<<"here";
            comment=1;
            counter=0;
            break;
        }else if(prev_residual>residual){
            prev_residual=residual;
            counter=0;
        }else{
            counter++;
        }
        if(counter>=10){
            comment=1;
            if(abs(wanted_over_relaxation-over_relaxation)<=0.01){
                line_gs=false;
            }
            // cout<<"here";
            break;
        }
        if(sym){
            iterations++;
            for(int i=temp.size()-2;i>=1;i--){
                coefficients=set_coeff(coefficients,temp,beta);
                rhs=set_rhs(rhs,temp,beta,i);

                apply_thomas_algorithm(coefficients,rhs);
                for(int j=0;j<temp[i].size();j++){
                    temp[i][j]+=(over_relaxation*(rhs[j]/coefficients[j][1]-temp[i][j]));
                }
            }
            residual=calculate_residual(temp,beta);
            if(abs(wanted_over_relaxation-over_relaxation)<=0.01){
                residual_vector_line.push_back(residual);
                // cout<<residual<<endl;
                if(prev_residual==residual || (prev_residual*1.1)<residual ){
                    line_gs=false;
                    comment=1;
                    cout<<"here"<<endl;
                    break;
                }
            }
            // cout<<residual<<" "<<iterations<<endl;
            
            if(prev_residual==residual || (prev_residual*1.1)<residual){
                comment=1;
                counter=0;
                // cout<<"here";
                break;
            }else if(prev_residual>residual){
                prev_residual=residual;
                counter=0;
            }else{
                counter++;
            }

            if(counter>=10){
            comment=1;
            if(abs(wanted_over_relaxation-over_relaxation)<=0.01){
                line_gs=false;
            }
            // cout<<"here";
            break;
            }
            // cout<<iterations<<" "<<residual<<endl;
        }
    }
    cout<<over_relaxation<<" "<<iterations<<" "<<comment<<endl;
    // cout<<residual_vector.size()<<endl;
}
void solve_using_ADI(vector<vector<float>>& temp,double beta,int max_iterations,double delta_x,double delta_y,vector<float>& residual_vector,int sym,float wanted_over_relaxation,bool & adi_converged,float over_relaxation){
    int iterations=0;
    int comment=0;
    double residual = calculate_residual(temp,beta); 
    double prev_residual=residual;
    int counter=0;
    if(abs(wanted_over_relaxation-over_relaxation)<=0.01){
        residual_vector.push_back(residual);
        // cout<<residual<<endl;
    }
    while(iterations<max_iterations && residual>=min((delta_x*delta_x),(delta_y*delta_y))){
        
        vector<vector<float>>coefficients(temp[0].size(),vector<float>(3,0));
        vector<float>rhs(temp[0].size(),0);
        // cout<<"hehe"<<endl;
        iterations++;
        for(int i=1;i<temp.size()-1;i++){
            coefficients=set_coeff(coefficients,temp,beta);
            rhs=set_rhs(rhs,temp,beta,i);
            apply_thomas_algorithm(coefficients,rhs);

            for(int j=0;j<temp[i].size();j++){
                temp[i][j]+=(over_relaxation*(rhs[j]/coefficients[j][1]-temp[i][j]));
            }
            // cout<<endl;
        }
        for(int i=1;i<temp[0].size()-1;i++){
            vector<vector<float>>coefficients_y(temp.size(),vector<float>(3,0));
            vector<float>rhs_y(temp.size(),0);
            
            coefficients_y=set_coeff_y(coefficients_y,temp,beta);
            rhs_y=set_rhs_y(rhs_y,temp,beta,i);
            apply_thomas_algorithm(coefficients_y,rhs_y);

            for(int j=0;j<temp.size();j++){
                temp[j][i]+=over_relaxation*(rhs_y[j]/coefficients_y[j][1]-temp[j][i]);
            }
            // cout<<endl;
        }
        
        residual=calculate_residual(temp,beta);
        if(abs(wanted_over_relaxation-over_relaxation)<=0.01){
            residual_vector.push_back(residual);
            if(prev_residual==residual || (prev_residual*1.1)<residual){
                // cout<<prev_residual<<" "<<residual<<" "<<over_relaxation<<endl;
                adi_converged=false;
                comment=1;
                break;
            }
        // cout<<residual<<endl;
        }
        // if(iterations>200){
        //     cout<<iterations<<" "<<prev_residual<<" "<<residual<<endl;
        // }
        if(prev_residual==residual || (prev_residual*1.1)<residual){
            comment=1;
            counter=0;
            break;
        }else if(prev_residual>residual){
            prev_residual=residual;
            counter=0;
        }else{
            counter++;
        }

        if(counter>=10){
            comment=1;
            if(abs(wanted_over_relaxation-over_relaxation)<=0.01){
                adi_converged=false;
            }
            break;
        }
        if(sym){
            iterations++;
        for(int i=temp.size()-2;i>=1;i--){
            coefficients=set_coeff(coefficients,temp,beta);
            rhs=set_rhs(rhs,temp,beta,i);
            apply_thomas_algorithm(coefficients,rhs);

            for(int j=0;j<temp[i].size();j++){
                temp[i][j]+=(over_relaxation*(rhs[j]/coefficients[j][1]-temp[i][j]));
            }
            // cout<<endl;
        }
        for(int i=temp[0].size()-2;i>=1;i--){
            vector<vector<float>>coefficients_y(temp.size(),vector<float>(3,0));
            vector<float>rhs_y(temp.size(),0);
            
            coefficients_y=set_coeff_y(coefficients_y,temp,beta);
            rhs_y=set_rhs_y(rhs_y,temp,beta,i);
            apply_thomas_algorithm(coefficients_y,rhs_y);

            for(int j=0;j<temp.size();j++){
                temp[j][i]+=over_relaxation*(rhs_y[j]/coefficients_y[j][1]-temp[j][i]);
            }
            // cout<<endl;
        }
        residual=calculate_residual(temp,beta);
        if(abs(wanted_over_relaxation-over_relaxation)<=0.01){
            residual_vector.push_back(residual);
            if(prev_residual==residual || (prev_residual*1.1)<residual){
                // cout<<prev_residual<<" "<<residual<<" "<<over_relaxation<<endl;
                adi_converged=false;
                comment=1;
                break;
            }
        // cout<<residual<<endl;
        }
        if(prev_residual==residual || (prev_residual*1.1)<residual){
            comment=1;
            counter=0;
            break;
        }else if(prev_residual>residual){
            prev_residual=residual;
            counter=0;
        }else{
            counter++;
        }
        if(counter>=10){
            comment=1;
            if(abs(wanted_over_relaxation-over_relaxation)<=0.01){
                adi_converged=false;
            }
            break;
        }
        }
        
    }
    cout<<over_relaxation<<" "<<iterations<<" "<<comment<<endl;
}

int main(){

    ifstream infile("input.dat");
    string key;
    float thermal_conductivity,thermal_diffusivity,min_over_relaxation=1,max_over_relaxation=1,step_over_relaxation=0.2,length =0.3,height =0.4,temp_x_0,temp_x_length,temp_y_0 ,temp_y_height,wanted_over_relaxation=1;
    double delta_x,delta_y;
    int max_iterations,symmetric=0;

    while (infile >> key) {
        if (key == "thermal_conductivity") {infile >> thermal_conductivity;
        }else if (key == "temp_x_0") {infile >> temp_x_0;
        }else if(key=="temp_x_length"){infile>>temp_x_length;
        }else if(key=="temp_y_0"){infile>>temp_y_0;
        }else if(key=="temp_y_height"){infile>>temp_y_height;
        }else if(key=="thermal_diffusivity"){infile>>thermal_diffusivity;
        }else if(key=="min_over_relaxation"){infile>>min_over_relaxation;
        }else if(key=="max_over_relaxation"){infile>>max_over_relaxation;
        }else if(key=="step_over_relaxation"){infile>>step_over_relaxation;
        }else if(key=="length"){infile>>length;
        }else if(key=="height"){infile>>height;
        }else if(key=="delta_x"){infile>>delta_x;
        }else if(key=="delta_y"){infile>>delta_y;
        }else if(key=="max_iter"){infile>>max_iterations;
        }else if(key=="wanted_over_relaxation"){infile>>wanted_over_relaxation;
        }else if(key=="symmetric"){infile>>symmetric;
        }
    }

    ll grid_points_y=(ll)(height/delta_y);
    ll grid_points_x=(ll)(length/delta_x);
    double beta=delta_x/delta_y;
    vector<int>comment(((int)((max_over_relaxation-min_over_relaxation+0.01)/step_over_relaxation)+1),0);

    // cout<<comment.size()<<endl;
    // cout<<(max_over_relaxation-min_over_relaxation+0.01)/step_over_relaxation +1<<endl;

    vector<vector<float>>test(grid_points_y+1,vector<float>(grid_points_x+1,0));
    int j=0;
    vector<float>residual_vector;
    int point_gs=1;
    cout<<"point GS"<<endl;
    for(float i=min_over_relaxation;i<=(max_over_relaxation+0.01);i+=step_over_relaxation){
        vector<vector<float>>temp(grid_points_y+1,vector<float>(grid_points_x+1,0));
        apply_boundary_conditions(temp,temp_x_0,temp_x_length,temp_y_0,temp_y_height);
        solve_using_point_GS(temp,i,beta,max_iterations,delta_x,delta_y,residual_vector,comment,j,wanted_over_relaxation,symmetric);
        
        if(i<=wanted_over_relaxation+0.01 && i>=wanted_over_relaxation-0.01){
            test=temp;
            if(comment[j]==1){
                point_gs=0;
            }
        }
        j++;
    }

    cout<<endl;

    bool line_gs=true;j=0;
    vector<float>residual_vector_line;
    cout<<"Line GS"<<endl;
    for(float i=min_over_relaxation;i<=(max_over_relaxation+0.01);i+=step_over_relaxation){
        vector<vector<float>>temp_line_gs(grid_points_y+1,vector<float>(grid_points_x+1,0));
        apply_boundary_conditions(temp_line_gs,temp_x_0,temp_x_length,temp_y_0,temp_y_height);
        solve_using_line_GS(temp_line_gs,beta,max_iterations,delta_x,delta_y,residual_vector_line,symmetric,wanted_over_relaxation,i,line_gs);
    }
    // cout<<residual_vector_line[255]<<" "<<residual_vector_line[254]<<endl;
    cout<<endl;
    bool adi_converged=true;j=0;
    vector<float>residual_vector_adi;
    cout<<"ADI"<<endl;
    for(float i=min_over_relaxation;i<=(max_over_relaxation+0.01);i+=step_over_relaxation){
        vector<vector<float>>temp_adi(grid_points_y+1,vector<float>(grid_points_x+1,0));
        apply_boundary_conditions(temp_adi,temp_x_0,temp_x_length,temp_y_0,temp_y_height);
        solve_using_ADI(temp_adi,beta,max_iterations,delta_x,delta_y,residual_vector_adi,symmetric,wanted_over_relaxation,adi_converged,i);
    }
    
    // for(int i=0;i<test.size();i++){
    //     for(int j=0;j<test[0].size();j++){
    //         cout<<test[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    cout<<endl;
    if(point_gs==1){
        cout<<"SOR"<< " : over_relaxation = "<<wanted_over_relaxation<<endl;
    }else{
        cout<<"SOR NOT CONVERGED "<< " : over_relaxation = "<<wanted_over_relaxation<<endl;
    }
    
    for(int i=0;i<residual_vector.size();i++){
        if(i%20==0){
            cout<<i<<" "<<residual_vector[i]<<endl;
        }
    }
    if((residual_vector.size()-1)%20 !=0){
        cout<<residual_vector.size()-1<<" "<<residual_vector[residual_vector.size()-1]<<endl;
    }

    residual_vector.clear();
    cout<<endl;

    if(line_gs){
        cout<<"Line GS CONVERGED"<<endl;
    }else{
        cout<<"Line GS NOT CONVERGED"<<endl;
    }
    for(int i=0;i<residual_vector_line.size();i++){
        if(i%20==0){
            cout<<i<<" "<<residual_vector_line[i]<<endl;
        }
    }
    if((residual_vector_line.size()-1)%20 !=0){
        cout<<residual_vector_line.size()-1<<" "<<residual_vector_line[residual_vector_line.size()-1]<<endl;
    }
    cout<<endl;
    if(adi_converged){
        cout<<"ADI CONVERGED"<<endl;
    }else{
        cout<<"ADI NOT CONVERGED"<<endl;
    }
    for(int i=0;i<residual_vector_adi.size();i++){
        if(i%20==0){
            cout<<i<<" "<<residual_vector_adi[i]<<endl;
        }
    }
    if((residual_vector_adi.size()-1)%20 !=0){
        cout<<residual_vector_adi.size()-1<<" "<<residual_vector_adi[residual_vector_adi.size()-1]<<endl;
    }







    // residual_vector.clear();
    // cout<<endl;

    // apply_boundary_conditions(temp,temp_x_0,temp_x_length,temp_y_0,temp_y_height);
    // residual_vector.clear();
    // bool adi_converged=true;
    // solve_using_ADI(temp,beta,max_iterations,delta_x,delta_y,residual_vector,symmetric,wanted_over_relaxation,adi_converged);
    // if(adi_converged){
    //     cout<<"ADI CONVERGED"<<endl;
    // }else{
    //     cout<<"ADI NOT CONVERGED"<<endl;
    // }
    
    // // cout<<residual_vector.size()<<endl;
    // for(int i=0;i<residual_vector.size();i++){
    //     if(i%20==0){
    //         cout<<i<<" "<<residual_vector[i]<<endl;
    //     }
    // }
    // if((residual_vector.size()-1)%20 !=0){
    //     cout<<residual_vector.size()-1<<" "<<residual_vector[residual_vector.size()-1]<<endl;
    // }

    return 0;
}