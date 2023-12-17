#include "main.h"

int numUpdates = 0; // current number of changes (add/remove Customers in min-heap restaurant)

class HuffNode{
public:
    char data;
    int weight;
    HuffNode *left, *right;

    HuffNode(char data = ' ', int weight = 0, HuffNode* left = nullptr, HuffNode* right = nullptr): data(data), weight(weight), left(left), right(right) {}
};

unordered_map<char, int> freqCounter(const string& input){
    unordered_map<char, int> result;
    for (char c: input) result[c]++;
    return result;
}
string caesarEncode(const string& name) {
    unordered_map<char, int> charFreq = freqCounter(name);
    string result = name;
    for (char &i: result) i = i / 32 * 32 + 1 + (i % 32 - 1 + charFreq[i]) % 26;
    return result;
}
int mod1024(string binStr) {
    string result;
    unsigned int len = binStr.length();
    for (unsigned int i = 1; i <= 10 && i <= len; i++) result += binStr[len - i];
    return stoi(result, nullptr, 2); // convert to base 2
}
void deleteHuffTree(HuffNode* root) {
    if (!root) return;
    deleteHuffTree(root->left);
    deleteHuffTree(root->right);
    delete root;
}
HuffNode* rotateLeft(HuffNode* root) {
    HuffNode *newRoot = root->right, *child = root->right->left;
    newRoot->left = root;
    root->right = child;
    return newRoot;
}
HuffNode* rotateRight(HuffNode* root) {
    HuffNode *newRoot = root->left, *child = root->left->right;
    newRoot->right = root;
    root->left = child;
    return newRoot;
}
int getHeight(HuffNode* root){ // a HuffTree has at most 103 nodes, it's okay to traverse
    if (!root) return 0;
    int leftH = getHeight(root->left), rightH = getHeight(root->right);
    return 1 + (leftH > rightH ? leftH : rightH);
}
int balanceFactor(HuffNode* root){
    if (!root) return 0;
    return getHeight(root->left)- getHeight(root->right);
}
HuffNode* balanceNode(HuffNode* root, bool &rotated) {
    if (!root || rotated) return root;
    if (balanceFactor(root) > 1) {
        if (balanceFactor(root->left) < 0) root->left = rotateLeft(root->left);
        root = rotateRight(root);
        rotated = true;
    } else if (balanceFactor(root) < -1) {
        if (balanceFactor(root->right) > 0) root->right = rotateRight(root->right);
        root = rotateLeft(root);
        rotated = true;
    }
    root->left = balanceNode(root->left, rotated);
    root->right = balanceNode(root->right, rotated);
    return root;
}
HuffNode* balanceTree(HuffNode* root) {
    for (int i = 0; i < 3; i++) {
        bool rotated = false;
        root = balanceNode(root, rotated);
    } // balance the tree 3 times, if it's still unbalanced, then it's okay to leave it
    return root;
}
string huffEncode(char c, HuffNode* root, string currString = ""){
    if (!root) return "";
    else if (root->data == c) return currString;
    else return !huffEncode(c, root->left, currString+"0").empty() ? huffEncode(c, root->left, currString+"0") : huffEncode(c, root->right, currString+"1");
}
void printHuffmanInOrder(HuffNode* root) {
    if (!root) return;
    printHuffmanInOrder(root->left);
    if (isalpha(root->data)) cout << root->data << endl; else cout << root->weight << endl;
    printHuffmanInOrder(root->right);
}
HuffNode* makeHuffTree(const string& input) {
    unordered_map<char, int> charFreq = freqCounter(input);
    vector<HuffNode *> huffQueue;
    for (auto i: charFreq) {
        auto *temp = new HuffNode(i.first, i.second);
        huffQueue.push_back(temp); // delete temp;
    }
    while (huffQueue.size() > 1) {
        stable_sort(begin(huffQueue), end(huffQueue), [](const auto &a, const auto &b) {
            return (a->weight < b->weight) || (a->weight == b->weight && (a->data / 32 > b->data / 32 || (a->data / 32 == b->data / 32 && a->data < b->data)));
        });
        HuffNode *leftT = huffQueue[0], *rightT = huffQueue[1];
        int totalWeight = leftT->weight + rightT->weight;
        huffQueue.erase(begin(huffQueue), begin(huffQueue) + 2);
        auto *temp = new HuffNode(' ', totalWeight, leftT, rightT);
        temp = balanceTree(temp);
        huffQueue.push_back(temp);
    }
    return huffQueue[0];
}

class Customer{
public:
    int result;
    Customer *left, *right;
    Customer(int res, Customer* l = nullptr, Customer* r = nullptr): result(res), left(l), right(r) {}
};

class Q { // same as Assignment 1, but left = prev, right = next for BST + LL practice
public:
    Customer *head, *tail;
    int count;

    Q() : head(nullptr), tail(nullptr), count(0) {}

    ~Q() {
        while (head) delete (dequeue());
        delete head;
    }

    void enqueue(Customer *c) { // add to tail
        if (!head) head = tail = c;
        else {
            c->left = tail;
            tail = tail->right = c;
        }
        count++;
    }

    Customer *dequeue() { // take the head of the queue, can be used to put in table
        Customer *pos = head;
        head = head->right;
        if (count == 1) tail = nullptr; else head->left = nullptr;
        pos->right = nullptr;
        count--;
        return pos;
    }
};

class BST {
public:
    Customer *root;
    Q *zoneTimeList;

    BST(Customer *root = nullptr) : root(root), zoneTimeList(new Q()) {}

    ~BST() {
        deleteHelper(root);
        delete zoneTimeList;
    }

    Customer *minValueNode(Customer *r) { return r->left ? minValueNode(r->left) : r; }

    Customer *addHelper(Customer *r, Customer *c) { // add customer c to tree r
        if (!r) return c;
        if (c->result < r->result) r->left = addHelper(r->left, c);
        else r->right = addHelper(r->right, c);
        return r;
    }

    void push(Customer *c) {
        root = addHelper(root, c);
        zoneTimeList->enqueue(new Customer(c->result));
    }

    void pop() { // remove front of the queue + remove from BST
        Customer *temp = zoneTimeList->dequeue();
        root = removeHelper(root, temp->result);
        delete temp;
    }

    Customer *removeHelper(Customer *r, int key) {
        if (!r) return r;
        if (key < r->result) r->left = removeHelper(r->left, key);
        else if (key > r->result) r->right = removeHelper(r->right, key);
        else if (!(r->left) || !(r->right)) { // 0 or 1 child
            Customer *temp = r->left ? r->left : r->right;
            delete r;
            return temp;
        } else {
            Customer *temp = minValueNode(r->right);
            r->result = temp->result;
            r->right = removeHelper(r->right, temp->result);
        }
        return r;
    }

    void deleteHelper(Customer *r) {
        if (!r) return;
        deleteHelper(r->left);
        deleteHelper(r->right);
        delete r;
    }
};

class GojoRes {
public:
    int numZones;
    BST **zones;

    GojoRes(int numZones) : numZones(numZones), zones(new BST *[numZones]) {
        for (int i = 0; i < numZones; i++)
            zones[i] = new BST();
    }

    ~GojoRes() {
        for (int i = 0; i < numZones; i++) delete zones[i];
        delete[]zones;
    }

    void add(Customer *c) const { zones[c->result % numZones]->push(c); }

    void printCustomersInOrder(Customer *r) { // LNR
        if (!r) return;
        printCustomersInOrder(r->left);
        cout << r->result << endl;
        printCustomersInOrder(r->right);
    }

    void printBSTInOrder(int idx) { printCustomersInOrder(zones[idx]->root); }

    int nCrMod(int n, int r) const {
        int nCr[n + 1][r + 1]; // dynamic programming
        for (int i = 0; i <= n; i++) {
            for (int j = 0; j <= r; j++)
                nCr[i][j] = (j == 0 || j == i ? 1 : (nCr[i - 1][j - 1] + nCr[i - 1][j])) % numZones;
        }
        return nCr[n][r];
    }

    void pushBSTPreOrder(Customer *root, vector<int> &v) {
        if (!root) return;
        v.push_back(root->result);
        pushBSTPreOrder(root->left, v);
        pushBSTPreOrder(root->right, v);
    }

    int numPermutationsBST(vector<int> &arr) {
        if (arr.size() <= 2) return 1 % numZones;
        vector<int> leftTree, rightTree;
        int root = arr[0];
        for (unsigned int i = 1; i < arr.size(); i++) {
            if (arr[i] < root) leftTree.push_back(arr[i]); else rightTree.push_back(arr[i]);
        } // modulo in each step
        return numPermutationsBST(leftTree) * numPermutationsBST(rightTree) % numZones *
               nCrMod(arr.size() - 1, leftTree.size()) % numZones;
    }

    void removeAllZones() {
        for (int i = 0; i < numZones; i++) {
            vector<int> v;
            pushBSTPreOrder(zones[i]->root, v);
            int numKicks = numPermutationsBST(v);
            for (int j = 0; j < numKicks && zones[i]->zoneTimeList->count; j++) zones[i]->pop();
        }
    }
};
class Compare {
public:
    bool operator()(tuple<int, int, int> a, tuple<int, int, int> b) { return a > b; }
};

class SukunaRes {
public:
    int numZones;
    Q **zones; // CAUTION: ZONE 0 = ID 1, +1 is just the extra, take the remainder is okay
    vector<tuple<int, int, int>> resHeap; // tuples of <numOfCustomers, lastUpdate[], zoneID>, use make_heap to heapify for this one

    SukunaRes(int numZones) : numZones(numZones), zones(new Q *[numZones]), resHeap({}) {
        for (int i = 0; i < numZones; i++)zones[i] = new Q();
    }

    ~SukunaRes() {
        for (int i = 0; i < numZones; i++) delete zones[i];
        delete[]zones;
        resHeap.clear();
    }

    int findInHeap(int resID) {
        tuple<int, int, int> res; // find the tuple of the correct zone
        for (unsigned int i = 0; i < resHeap.size(); i++) if (get<2>(resHeap[i]) == resID) return int(i);
        return -1;
    }

    void add(Customer *c) {
        int idx = findInHeap(c->result % numZones);
        if (idx < 0) resHeap.emplace_back(1, ++numUpdates, c->result % numZones);
        else resHeap[idx] = make_tuple(++get<0>(resHeap[idx]), ++numUpdates, get<2>(resHeap[idx]));
        make_heap(begin(resHeap), end(resHeap), [](const auto &a, const auto &b) { return a > b; });
        zones[c->result % numZones]->enqueue(c);
    }

    void removeMultipleZones(int num) {
        priority_queue<tuple<int, int, int>, vector<tuple<int, int, int>>, Compare> pq(begin(resHeap), end(resHeap));
        for (int i = 0; i < num && !pq.empty(); i++) {
            int resID = get<2>(pq.top());
            for (int j = 0; j < num && zones[resID]->count; j++) {
                Customer *temp = zones[resID]->dequeue();
                cout << temp->result << "-" << resID + 1 << endl;
                delete temp;
            }
            pq.pop();
            int idx = findInHeap(resID); // update the number of customers and last update time
            resHeap[idx] = make_tuple(zones[resID]->count, ++numUpdates, resID);
            make_heap(begin(resHeap), end(resHeap), [](const auto &a, const auto &b) { return a > b; });
            int idxFind = findInHeap(
                    resID); // after heapify, found the index of the zone again and delete if zone is empty
            if (get<0>(resHeap[idxFind]) == 0) {
                resHeap[idxFind].swap(resHeap[resHeap.size() - 1]);
                resHeap.pop_back();
                make_heap(begin(resHeap), end(resHeap), [](const auto &a, const auto &b) { return a > b; });
            }
        }
    }

    void printLIFO(int num, int idx) const { // print last num customers of zone idx (ID idx+1)
        Customer *temp = zones[idx]->tail;
        for (int i = 0; i < num && i < zones[idx]->count; i++, temp = temp->left)
            cout << idx + 1 << '-' << temp->result << endl;
    }

    void printMinHeapPreOrder(int num, int idx) { // left = 2*idx + 1, right = 2*idx+2
        if ((unsigned int) idx >= resHeap.size()) return;
        printLIFO(num, get<2>(resHeap[idx]));
        printMinHeapPreOrder(num, 2 * idx + 1);
        printMinHeapPreOrder(num, 2 * idx + 2);
    }
};

class TotalRes {
public:
    int size;
    HuffNode *lastHuffTree;
    GojoRes *resG;
    SukunaRes *resS;

    TotalRes(int maxSize) : size(maxSize), lastHuffTree(nullptr), resG(new GojoRes(maxSize)),
                            resS(new SukunaRes(maxSize)) {}

    ~TotalRes() {
        deleteHuffTree(lastHuffTree);
        delete resG;
        delete resS;
    }

    void LAPSE(const string &name) { // register name to a restaurant
        string res, encoded = caesarEncode(name);
        HuffNode *newHuff = makeHuffTree(encoded);
        if (freqCounter(name).size() < 3 || (isalpha(newHuff->data) && freqCounter(encoded).size() != 1)) {
            deleteHuffTree(newHuff);
            return;
        }
        deleteHuffTree(lastHuffTree);
        lastHuffTree = newHuff;
        for (unsigned int i = 0; i < name.length(); i++) res += huffEncode(encoded[i], newHuff);
        int result = res.empty() ? 0 : mod1024(res);
        auto *newC = new Customer(result);
        if (newC->result % 2 == 1) resG->add(newC); else resS->add(newC);
    }

    void KOKUSEN() const { resG->removeAllZones(); }

    void KEITEIKEN(int num) const { resS->removeMultipleZones(num); }

    void HAND() const { printHuffmanInOrder(lastHuffTree); }

    void LIMITLESS(int num) const { if (num > 0 && num <= size) resG->printBSTInOrder(num - 1); }

    void CLEAVE(int num) const { resS->printMinHeapPreOrder(num, 0); }
};

void simulate(const string& filename) {
    ifstream ss(filename);
    string str, maxsize, name, num;
    TotalRes *JJK;
    while (ss >> str) {
        if (str == "MAXSIZE") {
            ss >> maxsize;
            JJK = new TotalRes(stoi(maxsize));
        }
        else if (str == "LAPSE") {
            ss >> name;
            JJK->LAPSE(name);
        }
        else if (str == "KOKUSEN") { JJK->KOKUSEN(); }
        else if (str == "KEITEIKEN") {
            ss >> num;
            JJK->KEITEIKEN(stoi(num));
        }
        else if (str == "HAND") { JJK->HAND(); }
        else if (str == "LIMITLESS") {
            ss >> num;
            JJK->LIMITLESS(stoi(num));
        }
        else {
            ss >> num;
            JJK->CLEAVE(stoi(num));
        }
    }
    delete JJK;
}