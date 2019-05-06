//define a swap function to swap two variables.
void swap(int *xp, int *yp) 
{ 
    int temp = *xp; 
    *xp = *yp; 
    *yp = temp; 
} 
// T(n)=(n-1)*(n-1)*C
//An optimized version of Bubble Sort
void bubblesort(int arr[], int n) 
{ 
   int i, j; 
   bool swapped; 
   for (i = 1; i < n-1; i++) 
   { 
     //swapped = false; 
     for (j = 0; j < n-i-1; j++) 
     { 
        if (arr[j] > arr[j+1]) 
        { 
           swap(&arr[j], &arr[j+1]); 
           //swapped = true; 
        } 
     } 
  
     // IF no two elements were swapped by inner loop, then break 
     /*if (swapped == false) 
        break; */
     } 
} 

void InsertionSort(int arr[], int arr_size){
  if(arr_size > 1){
    int size = arr_size;
    //loop over all the elments except the first one
    for(int i = 1; i < size; ++i){
      int j = i - 1;
      int key = arr[i];
      /* Move elements of arrr[0, i-1], that are greater than key, 
      * to one position ahead of their current potion*/
      while(j >= 0 && arr[j] > key){
        //std::cout<<arr[j+1]<<"  "<<arr[j]<<std::endl;
        arr[j+1] = arr[j];
        //std::cout<<arr[j+1]<<"  "<<arr[j]<<std::endl;
        
        --j;
      }
      arr[j+1] = key;
    }//end of for loop
  }
}
 

//void Mergesort
//l -> left, m-> middle, r->right
// Merges two subarrays of arr[]. 
// // First subarray is arr[l..m] 
// // Second subarray is arr[m+1..r] 
void merge(int arr[], int l, int m, int r) 
{ 
    int i, j, k; //index for loop
    int n1 = m - l + 1; 
    int n2 =  r - m; 
    //length of the two sub arrays
    /* create temp arrays */
    int L[n1], R[n2]; 
    //copy the elements of arr into two temp arrays, L[] and R[]
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++) 
        L[i] = arr[l + i]; 
    for (j = 0; j < n2; j++) 
        R[j] = arr[m + 1+ j]; 

    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 

    while (i < n1 && j < n2) 
    { 
        if (L[i] <= R[j]) 
        { 
            arr[k] = L[i]; 
            i++; 
        } 
        else
        { 
            arr[k] = R[j]; 
            j++; 
        } 
        k++; 
    }   
   /* Copy the remaining elements of L[], if there 
 *        are any */
    while (i < n1) 
    { 
        arr[k] = L[i]; 
        i++; 
        k++; 
    } 
  
    /* Copy the remaining elements of R[], if there 
 *        are any */
    while (j < n2) 
    { 
        arr[k] = R[j]; 
        j++; 
        k++; 
    } 
} 


// call until l=r  
/* l is for left index and r is right index of the 
 *    sub-array of arr to be sorted */
void mergeSort(int arr[], int l, int r) 
{ 
    if (l < r) 
    { 
        // Same as (l+r)/2, but avoids overflow for 
        //         // large l and h 
       int m = l+(r-l)/2; 
       //calculate the middle index of the subarraies 
        // Sort first and second halves 
        mergeSort(arr, l, m); 
        mergeSort(arr, m+1, r); 
  
        merge(arr, l, m, r); 
    }        


}
//quickSort

/* Function to print an array */
void printArray(int arr[], int size) 
{ 
    int i; 
    for (i=0; i < size; i++) 
        printf("%d ", arr[i]); 
    //printf("n"); 
} 
//driver program to test the above functions
int sort_main() 
{ 
    int arr[] = {64, 34, 25, 12, 22, 11, 90};
    //get the total number of elements in an array 
    int n = sizeof(arr)/sizeof(arr[0]); 
    //bubblesort(arr, n); 
    //InsertionSort(arr, n);
    
    mergeSort(arr, 0, n-1);
    printf("Sorted array: \n"); 
    printArray(arr, n); 
    return 0; 
} 












