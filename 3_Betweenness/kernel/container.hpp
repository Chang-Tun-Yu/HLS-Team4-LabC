#include <iostream>
#include <iomanip>
#include <cstdlib>
using namespace std;

#define N 3534

typedef unsigned int data_t;
typedef bool flag_t;

// Class CFir definition
template<class data_T, class flag_T>
class Container {
 protected:
   data_T storage[N];
   data_T sp = 0; //start pointer
   data_T ep = 0; //end pointer
   data_T cnt = 0;

 private:
 public:
   flag_T   push_front(data_T element);
   flag_T   push_back(data_T element);

   data_T   pop();
   data_T   front();

   flag_T   empty();
   flag_T   full();

   data_T   start_pointer();
   data_T   end_pointer();
   data_T   get_element(data_T index);
};

template<class data_T, class flag_T>
flag_T Container<data_T, flag_T>::push_front(data_T element) {
    flag_T flag = ~full();
    if(flag)
    {
        storage[sp] = element;
        if(sp == 0)
            sp = N - 1;
        else
            sp = sp - 1;
        cnt = cnt + 1;
    }
    return flag;
}

template<class data_T, class flag_T>
flag_T Container<data_T, flag_T>::push_back(data_T element) {
    flag_T flag = ~full();
    if(flag)
    {
        storage[ep] = element;
        if(ep == N - 1)
            ep = 0;
        else
            ep = ep + 1;
        cnt = cnt + 1;
    }
    return flag;
}

template<class data_T, class flag_T>
data_T Container<data_T, flag_T>:: pop(){
    flag_T flag = ~empty();
    data_T data = -1;
    if(flag)
    {
        data = storage[ep];
        if(ep == 0)
            ep = N - 1;
        else
            ep = ep - 1;
        cnt = cnt - 1;
    }
    return data;
}

template<class data_T, class flag_T>
data_T Container<data_T, flag_T>:: front(){
    flag_T flag = ~empty();
    data_T data = -1;
    if(flag)
    {
        data = storage[sp];
        if(sp == N - 1)
            sp = 0;
        else
            sp = sp + 1;
        cnt = cnt - 1;
    }
    return data;
}

template<class data_T, class flag_T>
flag_T Container<data_T, flag_T>:: empty(){
    return (cnt == 0);
}

template<class data_T, class flag_T>
flag_T Container<data_T, flag_T>:: full(){
    return (cnt == N);
}

template<class data_T, class flag_T>
data_T Container<data_T, flag_T>:: start_pointer(){
    return sp;
}

template<class data_T, class flag_T>
data_T Container<data_T, flag_T>:: end_pointer(){
    return ep;
}

template<class data_T, class flag_T>
data_T Container<data_T, flag_T>:: get_element(data_T index){
    return storage[index];
}