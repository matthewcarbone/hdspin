#include <iostream>
#include <unordered_map> 

#ifndef LRU_H
#define LRU_H

// Modified from https://bhrigu.me/blog/2017/01/22/lru-cache-c-plus-plus-implementation/
// Uses a least-recently used queue
// Modified to use long long as keys and doubles as values

class Node
{
public:
    long long key;
    double value;
    Node *prev, *next;
    Node(long long k, double v): key(k), value(v), prev(NULL), next(NULL) {}
};


class DoublyLinkedList
{
private:
    Node *front, *rear;  
    bool isEmpty(){return rear == NULL;}

public:
    DoublyLinkedList(): front(NULL), rear(NULL) {}

    Node* add_page_to_head(long long key, double value)
    {
        Node *page = new Node(key, value);

        if (!front && !rear){front = rear = page;}
        else
        {
            page->next = front;
            front->prev = page;
            front = page;
        }

        return page;
    }

    void move_page_to_head(Node *page)
    {
        if (page == front) {return;}

        if (page == rear)
        {
            rear = rear->prev;
            rear->next = NULL;
        }
        else
        {
            page->prev->next = page->next;
            page->next->prev = page->prev;
        }

        page->next = front;
        page->prev = NULL;
        front->prev = page;
        front = page;
    }

    void remove_rear_page()
    {
        if(isEmpty()) {return;}
        if(front == rear)
        {
            delete rear;
            front = rear = NULL;
        }
        else
        {
            Node *temp = rear;
            rear = rear->prev;
            rear->next = NULL;
            delete temp;
        }
    }

    Node* get_rear_page() {return rear;}
  
};

class LRUCache
{
private:
    int capacity, size;
    DoublyLinkedList *pageList;
    std::unordered_map<long long, Node*> pageMap;

public:
    LRUCache(int capacity)
    {
        this->capacity = capacity;
        size = 0;
        pageList = new DoublyLinkedList();
        pageMap = std::unordered_map<long long, Node*>();
    }

    double get(long long key)
    {
        // If the key does not exist, return some number. It is extremely
        // unlikely that randomly, we will get exactly zero, so it is safe to
        // use it as the special value for a key not found.
        if(pageMap.find(key)==pageMap.end()) {return 0.0;}
        double val = pageMap[key]->value;

        // move the page to front
        pageList->move_page_to_head(pageMap[key]);
        return val;
    }

    void put(long long key, double value)
    {
        if(pageMap.find(key)!=pageMap.end())
        {
            // if key already present, update value and move page to head
            pageMap[key]->value = value;
            pageList->move_page_to_head(pageMap[key]);
            return;
        }

        if(size == capacity)
        {
          // remove rear page
          long long k = pageList->get_rear_page()->key;
          pageMap.erase(k);
          pageList->remove_rear_page();
          size--;
        }

        // add new page to head to Queue
        Node *page = pageList->add_page_to_head(key, value);
        size++;
        pageMap[key] = page;
    }

    ~LRUCache()
    {
        std::unordered_map<long long, Node*>::iterator i1;
        for(i1=pageMap.begin(); i1!=pageMap.end(); i1++)
        {
            delete i1->second;
        }
        delete pageList;
    }
};

#endif
