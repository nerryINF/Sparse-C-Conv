int num_sparse = 3; 

float features[3][5]={
{1,0,0,0,0},
{2,0,0,0,0},
{3,0,0,0,0}};

int indices[3][3]={
{1,1,1},
{3,3,4},
{6,3,6}};

int sparse_shape[3]={8,8,8};
int KL_VOL = 27; 

float kernel[27]={
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
1,
};