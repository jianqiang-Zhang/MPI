/*
 * 使用模拟退火算法(SA)求解TSP问题(以中国TSP问题为例)
 * 参考自《Matlab 智能算法30个案例分析》
 */
 #include"mpi.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>

#define T0 50000.0  // 初始温度
#define T_end (1e-8)
#define q  0.98   // 退火系数
#define L 1000  // 每个温度时的迭代次数，即链长
#define N 32  // 城市数量
//int city_list[N]; // 用于存放一个解

// 中国32个城市坐标
double city_pos[N][2] =
{
{1304,2312},{3639,1315},{4177,2244},{3712,1399},
{3488,1535},{3326,1556},{3238,1229},{4196,1004},
{4312,7901},{4386,5702},{3007,1970},{2562,1756},
{2788,1491},{2381,1676},{1332,6953},{3715,1678},
{3918,2179},{4061,2370},{3780,2212},{3676,2578},
{4029,2838},{4263,2931},{3429,1908},{3507,2367},
{3394,2643},{3439,3201},{2935,3240},{3140,3550},
{2545,2357},{2778,2826},{2370,2975},{4520,3412}

};

//函数声明
double distance(double *,double *); // 计算两个城市距离
double path_len(int* arr, double** city_divided, int CityNumber);  // 计算路径长度
void  init(int* city_list, int cityNumber);  //初始化函数
void create_new(int CityNumber, int* city_list); // 产生新解
// 距离函数
double distance(double * city1,double * city2)
{
    
    double x1 = *city1;
    double y1 = *(city1+1);
    double x2 = *(city2);
    double y2 = *(city2+1);
    double dis = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    return dis;
}

// 计算路径长度
double path_len(int * arr,double** city_divided,int CityNumber)
{

    //test
    if(CityNumber==N)
    {
	    //arr为某次结果的城市坐标序列
    //divided为数据坐标集
    double path = 0; // 初始化路径长度
    int index = *arr; // 定位到第一个数字(城市序号)
    
    
    for(int i=0;i< CityNumber -1;i++)
    {
        int index1 = *(arr+i);
        int index2 = *(arr+i+1);
        
        double dis = distance(city_pos[index1-1],
                                city_pos[index2-1]);
        path += dis;
    }
    int last_index = *(arr+ CityNumber -1); // 最后一个城市序号
    int first_index = *arr; // 第一个城市序号
    double last_dis = distance(city_pos[last_index-1],
                                city_pos[first_index-1]);
    path = path + last_dis;
    return path; // 返回总的路径长度
	   
    }
    else
    {
    //arr为某次结果的城市坐标序列
    //divided为数据坐标集
    double path = 0; // 初始化路径长度
    int index = *arr; // 定位到第一个数字(城市序号)
  
    for(int i=0;i< CityNumber -1;i++)
    {
        int index1 = *(arr+i);
        int index2 = *(arr+i+1);
        double dis = distance(city_divided[index1-1],
                                city_divided[index2-1]);
	
        path += dis;
    }
    int last_index = *(arr+ CityNumber -1); // 最后一个城市序号
    int first_index = *arr; // 第一个城市序号
    double last_dis = distance(city_divided[last_index-1],
                                city_divided[first_index-1]);
    path = path + last_dis;
    return path; // 返回总的路径长度
    }
}

// 初始化函数
void init(int *city_list,int cityNumber)
{
    for(int i=0;i< cityNumber;i++)
        city_list[i] = i+1;  // 初始化一个解
}

// 产生一个新解
// 此处采用随机交叉两个位置的方式产生新的解
void create_new(int CityNumber,int *city_list)
{
    double r1 = ((double)rand())/(RAND_MAX+1.0);
    double r2 = ((double)rand())/(RAND_MAX+1.0);
    int pos1 = (int)(CityNumber *r1); //第一个交叉点的位置
    int pos2 = (int)(CityNumber *r2);
    int temp = city_list[pos1];
    city_list[pos1] = city_list[pos2];
    city_list[pos2] = temp;   // 交换两个点
}

//用于全排列的两个函数
void swap(int i, int *a,int offset)
{
    int temp;
    temp = a[offset];
    a[offset] = a[i];
    a[i] = temp;
}

double perm(int* midlist,int** citylist, int* index, int len ,int cityNumber,int remain,int offset)
{
    //len为遍历次数
    //nn记录midlist的第一个index
    int n=0;//记录midList的第二个index
    static int nn=0;
    int temp=0;
    static double shortestListDistance;//记录最短的路径距离
    double distance;
    //int list[N];
    int* list=(int*)malloc(N*sizeof(int));
    if (offset == len - 1)
    {
    	if(remain==0)
    	{
    		for(int i=0;i<len;i++)
    		{
			for (int j = 0; j < cityNumber;j++ )
		        {
		            list[n] = citylist[index[i]][j];
		            //printf("nn=%d  , n=%d  List=%d\n",nn,n,list[n]);
		            n++;
		        }
		       	
    		}
    	}
    	else
    	{
		for (int i = 0; i < len; i++)
		{
		    if (index[i]<=(remain-1))     
		    {
		        for (int j = 0; j <cityNumber+1;j++ )
		        {
		            list[n] = citylist[index[i]][j];
		            n++;
		           
		        }
		    }
		    else
		    {
		        for (int j = 0; j < cityNumber; j++)
		        {
		        //printf("n=%d",n);
		            list[n] = citylist[index[i]][j];
		            n++;
		        }
		    }
		}
        }
        
        //判断是否是较短的路径
        if(nn==0)
        {
           shortestListDistance = path_len(list, city_pos, N);
           for(int i=0;i<N;i++)
           {
             midlist[i]=list[i];
           }
        }
        else
        {
           distance = path_len(list, city_pos, N);
           if(distance<=shortestListDistance)
           {
           	shortestListDistance=distance;
           	for(int i=0;i<N;i++)
           	{
             	  midlist[i]=list[i];
           	}
           }
           
        }
        nn++;
        return;
    }
    for (int i = offset; i < len; i++)
    {
        swap(i,index, offset);
        free(list);
        perm(midlist,citylist,index,len,cityNumber,remain,offset + 1);
        swap(i, index, offset);
    }
     return shortestListDistance;
}


// 主函数
int main(int argc, char* argv[])
{
    int myid, numprocs, namelen,cityNumber;//cityNumber是每个进程需要处理的城市数量
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    double** city_divided;
    double startwtime = 0.0,midwtime, endwtime;
    MPI_Status* status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Get_processor_name(processor_name, &namelen);

    if (numprocs > 1)
    {

        //分成相应的份数
        int remain = N % (numprocs - 1);
        if (remain != 0)
        {
            if (myid <= remain && myid != 0)
            {

                cityNumber = ((N - remain) / (numprocs - 1)) + 1;
                city_divided=(double**)malloc(cityNumber*sizeof(double));
                for (int i = 0; i < cityNumber; i++)
                {              
                    city_divided[i] = (double*)malloc(2*sizeof(double));
                }            
                for (int i = 0; i < cityNumber; i++)
                {
                    //printf("test,myid=%d\n",myid);
                    city_divided[i][0] = city_pos[cityNumber*(myid-1)+i][0];
                    city_divided[i][1] = city_pos[cityNumber*(myid-1)+i][1];
                }
            }
            else 
            {
                cityNumber = ((N - remain) / (numprocs - 1));          
                city_divided=(double**)malloc(cityNumber*sizeof(double));
                for (int i = 0; i < cityNumber; i++)
                {              
                    city_divided[i] = (double*)malloc(2*sizeof(double));
                }        
                for (int i = 0; i < cityNumber; i++)
                {
                    city_divided[i][0] = city_pos[cityNumber*(myid-1)+ remain+i][0];
                    city_divided[i][1] = city_pos[cityNumber*(myid-1)+ remain+i][1];
                }
            }

        }
        else
        {
            cityNumber = N / (numprocs - 1);
            city_divided=(double**)malloc(cityNumber*sizeof(double));
            for (int i = 0; i < cityNumber; i++)
            {              
               city_divided[i] = (double*)malloc(2*sizeof(double));
            }
            for (int i = 0; i < cityNumber; i++)
            {
                city_divided[i][0] = city_pos[cityNumber*(myid-1)+i][0];
                city_divided[i][1] = city_pos[cityNumber*(myid-1)+i][1];
            }
        }

        if (myid == 0)
        {
            int** cityList = (int**)malloc((numprocs-1) * sizeof(int));//用于存放子进程的结果
            
            if(remain==0)
            {
		for (int i = 0; i < numprocs-1; i++)
		{
	 		cityList[i] = (int*)malloc(cityNumber * sizeof(int));    	
		}
            }
            else
            {
            		//此时cityNumber为较小的数
		    for (int i = 0; i < numprocs-1; i++)
		    {
		    	
			if (i < remain)
			{
				cityList[i] = (int*)malloc((cityNumber+1) * sizeof(int));
			}
			else
			{
				cityList[i] = (int*)malloc(cityNumber * sizeof(int));
			}
		        
		    }
            }

            startwtime = MPI_Wtime();
            MPI_Barrier(MPI_COMM_WORLD);
            midwtime = MPI_Wtime();
            double cTime = midwtime - startwtime;//子节点的计算时间
            printf("slaver finished,time:%f s \t ,final linking...\n", cTime);
            //开始接收所有结果
            if(remain==0)
            {
            	    for (int i = 0; i < numprocs-1; i++)
		    {
		                      
		            MPI_Recv(cityList[i], cityNumber, MPI_INT, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }          
            }
            else
            {            
		    for (int i = 0; i < numprocs-1; i++)
		    {
		        if (i <= remain)
		        {
		            MPI_Recv(cityList[i], cityNumber+1, MPI_INT, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		        }
		        else
		        {                   
		            MPI_Recv(cityList[i], cityNumber, MPI_INT, i+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		        }
		    }
            }
            //开始将收到的list结果其连接        
            int index[numprocs - 1];//用于全排列,该数组的值作为citylist【n】的下标
            for (int i = 0; i < numprocs - 1; i++)
            {
                index[i] = i;
            }
            //int** midList;//全排列后的结果
            
            int allpai=numprocs-1;//全排列的解的数量
            for(int i=numprocs-2;i>0;i--)
            {
            	allpai=allpai*i;                       	
            }
  
            int bestList[N];//相对最优的路径
            double shortestListDistance1;
            shortestListDistance1 = perm(bestList, cityList, index, numprocs - 1, cityNumber, remain, 0);
            
            endwtime = MPI_Wtime();

            printf("\n\nall finished,shortest distance is :%f\t,total time :%f s\t\n",shortestListDistance1, (endwtime-startwtime));
            printf("relative best path :\n");
            for (int i = 0; i < N; i++)  // 输出最优路径
            {
                if(i==N-1)
                {
                	 printf("%d \n",bestList[i]);
                }
                else
                {
                	 printf("%d--->", bestList[i]);
                }
                
            }
        }
        else {
            
            int city_list[cityNumber];

            double T;
            int count = 0; // 记录降温次数
            T = T0; //初始温度
            init(city_list,cityNumber); //初始化一个解
            
             
            int city_list_copy[cityNumber]; // 用于保存原始解
            double f1, f2, df; //f1为初始解目标函数值，
                             //f2为新解目标函数值，df为二者差值
            double r; // 0-1之间的随机数，用来决定是否接受新解
            while (T > T_end) // 当温度低于结束温度时，退火结束
            {
                for (int i = 0; i < L; i++)
                {
                    // 复制数组
                                 
                    memcpy(city_list_copy, city_list, cityNumber * sizeof(int));
                    create_new(cityNumber,city_list); // 产生新解
                    f1 = path_len(city_list_copy, city_divided,cityNumber);
                    f2 = path_len(city_list, city_divided,cityNumber);
 
                    df = f2 - f1;//旧解减去新解时
                    // 以下是Metropolis准则
                    if (df >= 0)
                    {
                        //新生成的解更差时
                        r = ((double)rand()) / (RAND_MAX);
                        //printf("2222222\n");
                        if (exp(-df / T) <= r) // 保留原来的解
                        {
                            memcpy(city_list, city_list_copy, cityNumber * sizeof(int));
                        }
                    }
                }
                T *= q; // 降温
                count++;
            }
             
	            //将编号设置为全局的编号
            for (int i = 0; i < cityNumber; i++)
            {
                if (myid <= remain)
                {
                    //cityNumber较多的
                    city_list[i] += ((myid-1) * cityNumber);
                }
                else
                {
                    //cityNumber较少的
                    city_list[i] += ((myid-1) * cityNumber + remain);
                }
            }
            
            printf("No.%d proc which is on the mechine %s finished , the path is:\n",myid,processor_name);
            for (int i = 0; i < cityNumber; i++)
            {
                if (i == cityNumber - 1)
                {
                    printf("%d \n", city_list[i]);
                }
                else
                {
                    printf("%d -->", city_list[i]);
                }
            }


            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Send(city_list, cityNumber, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        int city_list[N];
        int count = 0; // 记录降温次数
        double T = T0; //初始温度
        init(city_list,N); //初始化一个解
        int city_list_copy[N]; // 用于保存原始解
        double f1, f2, df; //f1为初始解目标函数值，
                         //f2为新解目标函数值，df为二者差值
        double r; // 0-1之间的随机数，用来决定是否接受新解
        while (T > T_end) // 当温度低于结束温度时，退火结束
        {
            for (int i = 0; i < L; i++)
            {
                // 复制数组
                memcpy(city_list_copy, city_list, N * sizeof(int));
                create_new(N,city_list); // 产生新解
                f1 = path_len(city_list_copy,city_pos,N);
                f2 = path_len(city_list,city_pos,N);
                df = f2 - f1;
                // 以下是Metropolis准则
                if (df >= 0)
                {
                    //新生成的解更差时
                    r = ((double)rand()) / (RAND_MAX);
                    if (exp(-df / T) <= r) // 保留原来的解
                    {
                        memcpy(city_list, city_list_copy, N * sizeof(int));
                    }
                }
            }
            T *= q; // 降温
            count++;
        }
    }

    MPI_Finalize();
	
    return 0;
 }
    

