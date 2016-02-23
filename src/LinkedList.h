// Doubly linked. Assumes that T has members T *next and T *prev
template <typename T>
class LinkedList
{
public:

    LinkedList():_head{}, _tail{}
    {}

    void push_back(T *node)
    {
        if (!_head) {
            _head = _tail = node;
            node->prev = node->next = 0;
        } else {
            _tail->next = node;
            node->prev = _tail;
            node->next = 0;
        }
    }
    void insert_after(T *new_node, T *existing_node)
    {
        if (existing_node == _tail) {
            _tail = new_node;
            new_node->next = nullptr;
        } else {
            existing_node->next->prev = new_node;
            new_node->next = existing_node->next;
        }

        existing_node->next = new_node;
        new_node->prev = existing_node;
    }
    void insert_before(T *new_node, T *existing_node)
    {
        if (existing_node == _head)
        {
            _head = new_node;
            new_node->prev = nullptr;
        } else {
            existing_node->prev->next = new_node;
            new_node->prev = existing_node->prev;
        }
        existing_node->prev = new_node;
        new_node->next = existing_node;
    }


    T &head()
    {
        return *_head;
    }
    T &tail()
    {
        return *_tail;
    }
private:
    T *_head, *_tail;
};

// Doubly linked. Assumes that T has members T *next and T *prev
template <typename T>
class CircularList
{
public:

    CircularList() : _head {}
    {}

    void push_back(T *node)
    {
        if (!_head) {
            _head = node;
            node->prev = node->next = node;
        } else {
            node->prev = _head->prev;
            node->next = _head;
            _head->prev->next = node;
            _head->prev = node;

        }
    }
    void insert_after(T *new_node, T *existing_node)
    {
        new_node->prev = existing_node;
        new_node->next = existing_node->next;
        existing_node->next->prev = new_node;
        existing_node->next = new_node;
    }
    void insert_before(T *new_node, T *existing_node)
    {
        new_node->prev = existing_node->prev;
        new_node->next = existing_node;
        existing_node->prev->next = new_node;
        existing_node->prev = new_node;
    }

    T &head() {
        return *_head;
    }
    int length()
    {
        if (! _head ) return 0;
        int result = 0;
        T *node = _head;
        do {
            node = node->next;
            result ++;
        } while (node != _head);
        return result;
    }
    int numIntersects()
    {
        if (! _head ) return 0;
        int result = 0;
        T *node = _head;
        do {
            if (!node->intersect)
                result ++;
            node = node->next;
        } while (node != _head);
        return result;
    }

    ~CircularList()
    {
        T* element = _head;
        if (_head) {
            do {
                delete element;
            } while (element != _head);
        }
    }
private:
    T *_head;
};
