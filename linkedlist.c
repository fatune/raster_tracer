#include <stdlib.h>
/*
    initializing list

    NODE head;
    head.next = NULL;
    head.data = 0;
 
    another way to initialize list

    NODE * head;
    NODE headnode;
    head->next = NULL;
    InsertOrdered(head, 1);



    InsertOrdered(&head, i);

    Delete(&head, i);

void Traverse(NODE * head)
{
    NODE * current = head->next;
    while (current != NULL)
    {
        *******do stuff*********

        current = current->next;
    }
    printf("\n");
}

*/

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


void InsertOrdered(NODE *head, int data, long i, long j)
{
    NODE * newnode;
    newnode = (NODE *)malloc(sizeof(NODE));
    newnode->data = data;
    newnode->i = i;
    newnode->j = j;

    NODE * previous = head;
    NODE * current = head->next;
    while (current != 0 && data >= current->data)
    {
        previous = current;
        current = current->next;
    }
    previous->next = newnode;
 newnode->next = current;
}


int Delete(NODE *head, int data)
{
    NODE * previous = head;
    NODE * current = head->next;
    while (current != 0 && current->data != data)
    {
        previous = current;
        current = current->next;
    }
    if (current != head && current != 0) /* if list empty or data not found */
    {
        previous->next = current->next;
        free(current);
        return 0;
    }
    else
        return 1;
}

void Delete_second(NODE *head) {
    NODE * previous = head;
    NODE * current = head->next;
    previous->next = current->next;
    free(current);
}
