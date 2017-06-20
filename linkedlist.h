struct node {
    int data;
    struct node *next;
    long i;
    long j;
};

typedef struct node NODE;

void InsertOrdered(NODE *head, int data, long i, long j);
int Delete(NODE *head, int data);
void Delete_second(NODE *head);
