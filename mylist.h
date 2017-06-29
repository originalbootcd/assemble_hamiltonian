#ifndef MY_LIST_H
#define MY_LIST_H

template <class T> class List
{
public:
	List() {
		next = NULL;
		d = NULL;
	}

	List(List *list) {
		next = NULL;
		d = NULL;
		if (list->d)
			d = new T(*(list->d));
		if (list->next)
			next = new List(list->next);
	}

	List(List &list) {
		next = NULL;
		d = NULL;
		if (list.d)
			d = new T(*(list.d));
		if (list.next)
			next = new List(list.next);			
	}

	~List() {
		if (d) {
			delete d;
			d = NULL;
		}
		if (next) {
			delete next;
			next = NULL;
		}
	}
/*
	~List() {
		if (d) {
			free(d);
			d = NULL;
		}
		if (next) {
			free(next);
			next = NULL;
		}
	}
*/
	List<T> *next;
	T *d;
};
 ////!!!!!!  OPTIMIZED !!!!!!!!!!
template <class T> void Add(T *data, List<T> **head) {
	List<T> *newElement = new List<T>;
	//List<T> *newElement = (List<T>*)malloc(sizeof(List<T>));
	newElement->d = data;												//	newElement->d = new T(*data);
	newElement->next = *head;
	*head = newElement;
}


//////////// OPTIMIZED right order ////////////////
template <class T> void Add_tail(T *data, List<T> **head, List<T> **tail) {
	List<T> *newElement = new List<T>;
	newElement->d = data;	                                            //	newElement->d = new T(*data);
	newElement->next = NULL;

	if (*head) {
		(*tail)->next = newElement;
		*tail = (*tail)->next;
	}
	else {
		*head = newElement;
		*tail = newElement;
	}
}

/* //////// NONOPTIMIZED OLD ////////////////////
template <class T> void Add(T *data, List<T> **head) {
	List<T> *newElement = new List<T>;
	newElement->d = data;	                                            //	newElement->d = new T(*data);
	newElement->next = NULL;

	if (*head) {
		List<T> * tail = *head;
		while (tail->next)
			tail = tail->next;
		tail->next = newElement;
	}
	else
		*head = newElement;
}
*/

// wrong order (tail-head)
template <class T> void Del(List<T> **head) {
	if (*head) {
		List<T> *itemToBeDel = *head;
		*head = (*head)->next;
		itemToBeDel->next = NULL;
        itemToBeDel->d = NULL;
		delete itemToBeDel;
	}
	else cout << "Del() error -- the List is already empty" << endl; // error handling
}


#endif
