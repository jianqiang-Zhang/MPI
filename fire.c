/*
 * ʹ��ģ���˻��㷨(SA)���TSP����(���й�TSP����Ϊ��)
 * �ο��ԡ�Matlab �����㷨30������������
 */
 #include"mpi.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>

#define T0 50000.0  // ��ʼ�¶�
#define T_end (1e-8)
#define q  0.98   // �˻�ϵ��
#define L 1000  // ÿ���¶�ʱ�ĵ���������������
#define N 32  // ��������
//int city_list[N]; // ���ڴ��һ����

// �й�32����������
double city_pos[N][2] =
{
{1304,2312},{3639,1315},{4177,2244},{3712,1399},
{3488,1535},{3326,1556},{3238,1229},{4196,1004},
{4312,7901},{4386,5702},{3007,1970},{2562,1756},
{2788,1491},{2381,1676},{1332,695},
{3715,1678},{3918,2179},{4061,2370},
{3780,2212},{3676,2578},{4029,2838},
{4263,2931},{3429,1908},{3507,2367},
{3394,2643},{3439,3201},{2935,3240},
{3140,3550},{2545,2357},{2778,2826},
{2370,2975},{2323,2233} 
};

//��������
double distance(double *,double *); // �����������о���
double path_len(int* arr, double** city_divided, int CityNumber);  // ����·������
void  init(int* city_list, int cityNumber);  //��ʼ������
void create_new(int CityNumber, int* city_list); // �����½�
// ���뺯��
double distance(double * city1,double * city2)
{
    
    double x1 = *city1;
    double y1 = *(city1+1);
    double x2 = *(city2);
    double y2 = *(city2+1);
    double dis = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    return dis;
}

// ����·������
double path_len(int * arr,double** city_divided,int CityNumber)
{

    //test
    if(CityNumber==N)
    {
	    //arrΪĳ�ν���ĳ�����������
    //dividedΪ�������꼯
    double path = 0; // ��ʼ��·������
    int index = *arr; // ��λ����һ������(�������)
    
    
    for(int i=0;i< CityNumber -1;i++)
    {
        int index1 = *(arr+i);
        int index2 = *(arr+i+1);
        
        double dis = distance(city_pos[index1-1],
                                city_pos[index2-1]);
        path += dis;
    }
    int last_index = *(arr+ CityNumber -1); // ���һ���������
    int first_index = *arr; // ��һ���������
    double last_dis = distance(city_pos[last_index-1],
                                city_pos[first_index-1]);
    path = path + last_dis;
    return path; // �����ܵ�·������
	   
    }
    else
    {
    //arrΪĳ�ν���ĳ�����������
    //dividedΪ�������꼯
    double path = 0; // ��ʼ��·������
    int index = *arr; // ��λ����һ������(�������)
  
    for(int i=0;i< CityNumber -1;i++)
    {
        int index1 = *(arr+i);
        int index2 = *(arr+i+1);
        double dis = distance(city_divided[index1-1],
                                city_divided[index2-1]);
	
        path += dis;
    }
    int last_index = *(arr+ CityNumber -1); // ���һ���������
    int first_index = *arr; // ��һ���������
    double last_dis = distance(city_divided[last_index-1],
                                city_divided[first_index-1]);
    path = path + last_dis;
    return path; // �����ܵ�·������
    }
}

// ��ʼ������
void init(int *city_list,int cityNumber)
{
    for(int i=0;i< cityNumber;i++)
        city_list[i] = i+1;  // ��ʼ��һ����
}

// ����һ���½�
// �˴����������������λ�õķ�ʽ�����µĽ�
void create_new(int CityNumber,int *city_list)
{
    double r1 = ((double)rand())/(RAND_MAX+1.0);
    double r2 = ((double)rand())/(RAND_MAX+1.0);
    int pos1 = (int)(CityNumber *r1); //��һ��������λ��
    int pos2 = (int)(CityNumber *r2);
    int temp = city_list[pos1];
    city_list[pos1] = city_list[pos2];
    city_list[pos2] = temp;   // ����������
}

//����ȫ���е���������
void swap(int i, int *a,int offset)
{
    int temp;
    temp = a[offset];
    a[offset] = a[i];
    a[i] = temp;
}

double perm(int* midlist,int** citylist, int* index, int len ,int cityNumber,int remain,int offset)
{
    //lenΪ��������
    //nn��¼midlist�ĵ�һ��index
    int n=0;//��¼midList�ĵڶ���index
    static int nn=0;
    int temp=0;
    static double shortestListDistance;//��¼��̵�·������
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
        
        //�ж��Ƿ��ǽ϶̵�·��
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


// ������
int main(int argc, char* argv[])
{
    int myid, numprocs, namelen,cityNumber;//cityNumber��ÿ��������Ҫ����ĳ�������
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

        //�ֳ���Ӧ�ķ���
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
            int** cityList = (int**)malloc((numprocs-1) * sizeof(int));//���ڴ���ӽ��̵Ľ��
            
            if(remain==0)
            {
		for (int i = 0; i < numprocs-1; i++)
		{
	 		cityList[i] = (int*)malloc(cityNumber * sizeof(int));    	
		}
            }
            else
            {
            		//��ʱcityNumberΪ��С����
		    for (int i = 0; i < numprocs-1; i++)
		    {
		    	
			if (i <= remain)
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
            double cTime = midwtime - startwtime;//�ӽڵ�ļ���ʱ��
            printf("slaver finished,time:%fms \t ,final linking...\n", cTime * 1000);
            //��ʼ�������н��
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
            //��ʼ���յ���list���������        
            int index[numprocs - 1];//����ȫ����,�������ֵ��Ϊcitylist��n�����±�
            for (int i = 0; i < numprocs - 1; i++)
            {
                index[i] = i;
            }
            //int** midList;//ȫ���к�Ľ��
            
            int allpai=numprocs-1;//ȫ���еĽ������
            for(int i=numprocs-2;i>0;i--)
            {
            	allpai=allpai*i;                       	
            }
  
            int bestList[N];//������ŵ�·��
            double shortestListDistance1;
            shortestListDistance1 = perm(bestList, cityList, index, numprocs - 1, cityNumber, remain, 0);
            
            endwtime = MPI_Wtime();

            printf("all finished,shortest distance is :%f\t,total time :%fms\t\n",shortestListDistance1, (endwtime-startwtime) * 1000);
            printf("relative best path :\n");
            for (int i = 0; i < N; i++)  // �������·��
            {
                if(i==N-1)
                {
                	 printf("%d",bestList[i]);
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
            int count = 0; // ��¼���´���
            T = T0; //��ʼ�¶�
            init(city_list,cityNumber); //��ʼ��һ����
            
             
            int city_list_copy[cityNumber]; // ���ڱ���ԭʼ��
            double f1, f2, df; //f1Ϊ��ʼ��Ŀ�꺯��ֵ��
                             //f2Ϊ�½�Ŀ�꺯��ֵ��dfΪ���߲�ֵ
            double r; // 0-1֮�������������������Ƿ�����½�
            while (T > T_end) // ���¶ȵ��ڽ����¶�ʱ���˻����
            {
                for (int i = 0; i < L; i++)
                {
                    // ��������
                                 
                    memcpy(city_list_copy, city_list, cityNumber * sizeof(int));
                    create_new(cityNumber,city_list); // �����½�
                    f1 = path_len(city_list_copy, city_divided,cityNumber);
                    f2 = path_len(city_list, city_divided,cityNumber);
 
                    df = f2 - f1;//�ɽ��ȥ�½�ʱ
                    // ������Metropolis׼��
                    if (df >= 0)
                    {
                        //�����ɵĽ����ʱ
                        r = ((double)rand()) / (RAND_MAX);
                        //printf("2222222\n");
                        if (exp(-df / T) <= r) // ����ԭ���Ľ�
                        {
                            memcpy(city_list, city_list_copy, cityNumber * sizeof(int));
                        }
                    }
                }
                T *= q; // ����
                count++;
            }
             
	            //���������Ϊȫ�ֵı��
            for (int i = 0; i < cityNumber; i++)
            {
                if (myid <= remain)
                {
                    //cityNumber�϶��
                    city_list[i] += ((myid-1) * cityNumber);
                }
                else
                {
                    //cityNumber���ٵ�
                    city_list[i] += ((myid-1) * cityNumber + remain);
                }
            }
            
            printf("No  %d proc which is on the mechine %s finished , the path is:\n",myid,processor_name);
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
        int count = 0; // ��¼���´���
        double T = T0; //��ʼ�¶�
        init(city_list,N); //��ʼ��һ����
        int city_list_copy[N]; // ���ڱ���ԭʼ��
        double f1, f2, df; //f1Ϊ��ʼ��Ŀ�꺯��ֵ��
                         //f2Ϊ�½�Ŀ�꺯��ֵ��dfΪ���߲�ֵ
        double r; // 0-1֮�������������������Ƿ�����½�
        while (T > T_end) // ���¶ȵ��ڽ����¶�ʱ���˻����
        {
            for (int i = 0; i < L; i++)
            {
                // ��������
                memcpy(city_list_copy, city_list, N * sizeof(int));
                create_new(N,city_list); // �����½�
                f1 = path_len(city_list_copy,city_pos,N);
                f2 = path_len(city_list,city_pos,N);
                df = f2 - f1;
                // ������Metropolis׼��
                if (df >= 0)
                {
                    //�����ɵĽ����ʱ
                    r = ((double)rand()) / (RAND_MAX);
                    if (exp(-df / T) <= r) // ����ԭ���Ľ�
                    {
                        memcpy(city_list, city_list_copy, N * sizeof(int));
                    }
                }
            }
            T *= q; // ����
            count++;
        }
    }

    MPI_Finalize();
	
    return 0;
 }
    
