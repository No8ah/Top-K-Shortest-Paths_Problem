#include <iostream>
#include <iomanip>
#include<time.h>
#include <vector>
// 为了存储节点 --> 拿一个集合库来用
#include <unordered_set>
// 为了Ban掉节点的边的权值的入/出顶点 --> 拿一个Hash表进行使用
#include <unordered_map>
using namespace std;

#define MAX_DISTANCE  10000

template <class T>
class VecList{
    private:
        int capacity;
        int length;
        T* arr;
        void doubleListSize(){
            T * oldArr = arr;
            arr = new T[2*capacity];
            capacity = 2 * capacity;
            for(int i=0;i<length;i++){
                arr[i] = oldArr[i];
            }
            delete [] oldArr;
        }
    public:
        VecList(){
            length = 0;
            capacity = 100;
            arr = new T[capacity];
        }
        VecList(T* a, int n){
            length = n;
            capacity = 100 + 2*n;
            arr = new T[capacity];
            for(int i=0;i<n;i++){
                arr[i] = a[i];
            }
        }
        ~VecList(){
            delete [] arr;
        }
        int getLength(){
            return length;
        }
        bool isEmpty(){
            return length==0;
        }
        void insertEleAtPos(int i, T x){
            if(length==capacity)
                doubleListSize();
            if(i > length || i < 0)
                throw "Illegal position";
            for(int j=length;j>i;j--)
                arr[j] = arr[j-1];
            arr[i] = x;
            length++;
        }
        T deleteEleAtPos(int i){
            if(i >= length || i < 0)
                throw "Illegal position";
            T tmp = arr[i];
            for(int j=i;j<length-1;j++)
                arr[j] = arr[j+1];
            length--;
            return tmp;
        }
        void setEleAtPos(int i, T x){
            if(i >= length || i < 0)
                throw "Illegal position";
            arr[i] = x;
        }
        T getEleAtPos(int i){
            if(i >= length || i < 0)
                throw "Illegal position";
            return arr[i];
        }
        int locateEle(T x){
            for(int i=0;i<length;i++){
                if(arr[i]==x)
                    return i;
            }
            return -1;
        }
        void insertLast(T x){
            insertEleAtPos(length,x);
        }
        void printList(){
            for(int i=0;i<length;i++)
                cout << arr[i] << " ";
        }
        void reverseEleInVecList(){
            int front = 0;
            int end = length-1;
            while(front < end){
                T tmp = arr[front];
                arr[front] = arr[end];
                arr[end] = tmp;
                front++;
                end--;
            }
        }
};

template <class T>
struct DNode{
    T data;
    DNode<T>* next;
};

template <class T>
class LinkStack{
    private:
        DNode<T>* top;
        int length;
    public:
        LinkStack(){
            top = NULL;
            length = 0;
        }
        ~LinkStack(){
            while(top!=NULL){
                DNode<T>* tmp = top;
                top = top->next;
                delete tmp;
            }
        }
        int getLength(){
            return length;
        }
        bool isEmpty(){
            return length==0;
        }
        void push(T x){
            DNode<T>* tmp = new DNode<T>;
            tmp->data = x;
            tmp->next = top;
            top = tmp;
            length++;
        }
        T pop(){
            if(length==0) throw "Stack Empty!";
            DNode<T>* tmp = top;
            top = top->next;
            T tmpData = tmp->data;
            delete tmp;
            length--;
            return tmpData;
        }
        T getTop(){
            if(length==0) throw "Stack Empty!";
            return top->data;
        }
        void printStack(){
            cout << "Stack top: ";
            DNode<T>* tmp = top;
            while(tmp!=NULL){
                cout << tmp->data << " ";
                tmp = tmp->next;
            }
            cout << ":stack bottom" << endl;
        }
};

template <class T>
class LinkQueue{
    private:
        DNode<T>* front;
        DNode<T>* back;
        int length;
    public:
        LinkQueue(){
            front = new DNode<T>;
            front->next = NULL;
            back = front;
            length = 0;
        }
        ~LinkQueue(){
            while(front!=NULL){
                back = front;
                front = front->next;
                delete back;
            }
        }
        int getLength(){
            return length;
        }
        bool isEmpty(){
            return length==0;
        }
        void enQueue(T x){
            DNode<T>* tmpN = new DNode<T>;
            tmpN->data = x;
            tmpN->next = NULL;
            back->next = tmpN;
            back = tmpN;
            length++;
        }
        T deQueue(){
            if(length==0) throw "Queue Empty!";
            DNode<T>* tmpN = front->next;
            front->next = tmpN->next;
            T tmpD = tmpN->data;
            delete tmpN;
            length--;
            if(length==0) back = front;
            return tmpD;
        }
        T peekQueue(){
            if(length==0) throw "Queue Empty!";
            return front->next->data;
        }
        void printQueue(){
            cout << "Front of queue: ";
            DNode<T>* tmp = front->next;
            while(tmp!=NULL){
                cout << tmp->data << " ";
                tmp = tmp->next;
            }
            cout << ": back of queue" << endl;
        }
};

template <class T>
struct Edge{
    T start;
    T end;
    int weight;
};


template <class T>
class AMGraph{ 
    private:
        int numVer, numEdge;
        VecList<T> verList; 
        int** adjMatrix;
        bool directed;
        // 一个集合 --> used to collect vertices that are put into the black list
        unordered_set<T> black_list;
        unordered_map<T, unordered_map<T, int>> original_banned_vertex_edges;
        // 设置一下最大中继点个数
        int max_transition_vertex_num;

        void BFShelper(int st, bool* visited){
        	visited[st] = true;
        	cout << verList.getEleAtPos(st) << " ";
        	LinkQueue<int> q;
        	q.enQueue(st);
        	while(!q.isEmpty()){
        		int tmp = q.deQueue();
        		for(int k=0;k<numVer;k++){
        			// investigate adjMatrix[tmp][k]
        			if(adjMatrix[tmp][k]==0) continue;
        			if(visited[k]) continue;
        			visited[k] = true;
        			cout << verList.getEleAtPos(k) << " ";
        			q.enQueue(k);
				}
			}
		}

		void DFShelper(int st, bool* visited){
			visited[st] = true;
			cout << verList.getEleAtPos(st) << " ";
			for(int k=0;k<numVer;k++){
				if(adjMatrix[st][k]==0) continue;
				if(visited[k]) continue;
				DFShelper(k,visited);
			}
		}
    public:
        AMGraph(){ // we don't want this used.
        }
        AMGraph(T* arr, int n, bool dir=false){
            // n for number of vertices
            // default for undirected graph
            // edges to be added later.

            numVer = n;
			numEdge = 0;
            max_transition_vertex_num = n;
			for(int i=0;i<n;i++){
				verList.insertLast(arr[i]);
			}
			directed = dir;

			adjMatrix = new int*[n];
			for(int i=0;i<n;i++){
				adjMatrix[i] = new int[n];
				for(int j=0;j<n;j++)
					adjMatrix[i][j] = 0;
			}
        }
        AMGraph(T* arr, int n, Edge<T>** eArr, int e,bool dir=false){
            // n for number of vertices
            // default for undirected graph
            // edges to be added now.

            numVer = n;
			numEdge = 0;
			for(int i=0;i<n;i++){
				verList.insertLast(arr[i]);
			}
			directed = dir;

			adjMatrix = new int*[n];
			for(int i=0;i<n;i++){
				adjMatrix[i] = new int[n];
				for(int j=0;j<n;j++)
					adjMatrix[i][j] = 0;
			}

			for(int i=0;i<e;i++){
				addEdge(eArr[i]->start,eArr[i]->end, eArr[i]->weight);
			}
        }
        ~AMGraph(){
            for(int i=0;i<numVer;i++){
            	delete [] adjMatrix[i];
			}
			delete [] adjMatrix;
        }
        void addEdge(Edge<T> e){
            addEdge(e.start,e.end, e.weight);
        }
        void addEdge(T st, T en, int weight){
            int sIndex = verList.locateEle(st);
            int eIndex = verList.locateEle(en);
            if(adjMatrix[sIndex][eIndex]!=0) return;
            numEdge++;
            adjMatrix[sIndex][eIndex] = weight;
            if(!directed) adjMatrix[eIndex][sIndex] = weight;
        }
        void removeEdge(Edge<T> e){
            removeEdge(e.start,e.end);
        }
        void removeEdge(T st, T en){
            int sIndex = verList.locateEle(st);
            int eIndex = verList.locateEle(en);
            if(adjMatrix[sIndex][eIndex]==0) return;
            numEdge--;
            adjMatrix[sIndex][eIndex] = 0;
            if(!directed) adjMatrix[eIndex][sIndex] = 0;
        }
        void printGraph(){
            cout << "Vertices:" << endl;
            for(int i=0;i<numVer;i++)
                cout << verList.getEleAtPos(i) << " ";

            cout << endl << "Edges:" << endl;
            char sLeft = (directed ? '<' : '(');
            char sRight = (directed ? '>' : ')');
            for(int i=0;i<numVer;i++){
                for(int j=i+1;j<numVer;j++){
                    if(adjMatrix[i][j] != 0) cout << sLeft << verList.getEleAtPos(i) << "," << verList.getEleAtPos(j) << sRight << ", weight=" << adjMatrix[i][j] <<endl;
                }
            }
            if(!directed) return;
            for(int i=0;i<numVer;i++){
                for(int j=0;j<i;j++){
                    if(adjMatrix[i][j] != 0) cout << sLeft << verList.getEleAtPos(i) << "," << verList.getEleAtPos(j) << sRight << ", weight=" << adjMatrix[i][j] <<endl;
                }
            }
        }
        int** getMatrix(){
            return adjMatrix;
        }
        //---------------------------------------------------------- 核心功能函数 ----------------------------------------------------------//

            // ---------------------------------------- Ban ---------------------------------------- //
                void addVertexIntoBlackList(T Ban_Vertex){
                    // 此时顶点不在图里面 --> throw error and return 
                    int ban_vertex_index = verList.locateEle(Ban_Vertex);
                    if(ban_vertex_index==-1){
                        cout << "Vertex " << Ban_Vertex << " is not in Graph, please input another vertex again." << endl;
                        return;
                    }
                    // 此时顶点在图里面 --> 将顶点丢进黑名单里面 并将它的入边/顶点和出边/顶点 记录保存 并 禁用
                    else{
                        cout << "Vertex " << Ban_Vertex << " is added to black list" << endl;
                        // 将顶点放入到黑名单里面
                        black_list.insert(Ban_Vertex);
                        // 将禁用顶点的入/出顶点 和 边的权值 保存下来
                        int ban_vertex_index = verList.locateEle(Ban_Vertex);
                        saveInVertexAndEdgeWeights(ban_vertex_index,Ban_Vertex);
                        saveOutVertexAndEdgeWeights(ban_vertex_index,Ban_Vertex);
                    }
                }
            // ---------------------------------------- Ban ---------------------------------------- //

            // ----------------------------------------  Unban ---------------------------------------- //
                void removeVertexOutOfBlackList(T Unban_Vertex){
                    int unban_vertex_index = verList.locateEle(Unban_Vertex);
                    // 此时顶点不在黑名单里面 或者 顶点不在图里面 --> throw error and return
                    if(unban_vertex_index == -1){
                        cout << "Vertex " << Unban_Vertex << " is not in graph, please input other command again." << endl;
                        return;
                    }
                    if(!isVertexInBlackList(Unban_Vertex)){
                        cout << "Vertex " << Unban_Vertex << " is not in black list, please input other command again." << endl;
                        return;
                    }
                    // 此时顶点既在图里面 又 在黑名单里面 --> 将顶点从黑名单中移除 并将它的入边/顶点和出边/顶点 恢复
                    else{
                        // 将禁用顶点的入边/顶点 和 出边/顶点 恢复
                        recoverVertexAndEdgeWeight(Unban_Vertex, unban_vertex_index);
                        // 将禁用顶点移出黑名单
                        black_list.erase(Unban_Vertex);
                    }
                }
            // ----------------------------------------  Unban ---------------------------------------- //

            // ----------------------------------------  MaxTrans ---------------------------------------- //
                void limitMaxTransitVertexAs(int Max_Trans_Num){
                    //   0 大于等于 最大中继点的个数 小于等于 numVer-2 
                    if(Max_Trans_Num >= 0 && Max_Trans_Num <= numVer-2){
                        cout << "Max transition vertex numver has been set: " << Max_Trans_Num << endl;
                        max_transition_vertex_num = Max_Trans_Num;
                    }
                    else if(Max_Trans_Num == -1){
                        cout << "Max transition vertex number is now unlimited." << endl;
                        max_transition_vertex_num = numVer;
                    }
                    else{
                        cout << "Please input once again." << endl;
                    }
                }
            // ----------------------------------------  MaxTrans ---------------------------------------- //

            // ----------------------------------------  Paths ---------------------------------------- //
                void findTopKShortestPath(T Start_Vertex, T End_Vertex, int Top_K_Shortest_Path_Num){
                    // 初始化
                    int* shortest_distance = new int[numVer]; 
                    int* pre_vertex_index = new int[numVer];
                    bool* visited = initVisitedArray();
                    bool* is_in_set = new bool[numVer];
                    int* transit_vertex_counter = new int[numVer];
                    initDijkstra(shortest_distance,pre_vertex_index,transit_vertex_counter,is_in_set);
                    // 开始顶点到自身的距离当然设置为0
                    // Dijkstra : 首先算出最短路径
                    Dijkstra(Start_Vertex,shortest_distance,pre_vertex_index,is_in_set,transit_vertex_counter);
                    // Dijkstra : 首先算出最短路径
                    // 用于存储K条最短路径的向量表 : store go through vertices & store path weight
                    VecList<PathInformation*>* top_k_shortest_path = new VecList<PathInformation*>;
                    // 用于存储K条最短路径的向量表 : store go through vertices & store path weight

                    // 首先找出第一条最短路径
                    VecList<T>* first_shortest_path = findFirstShortestPath(Start_Vertex,End_Vertex,shortest_distance,pre_vertex_index);
                    // 首先找出第一条最短路径
                    // 如果可以找到第一条最短路径 --> 放入top k shortest path里面
                    int first_shortest_path_vertex_num = first_shortest_path->getLength();
                    int End_Vertex_Index = verList.locateEle(End_Vertex);
                    if(first_shortest_path_vertex_num > 0
                    && shortest_distance[End_Vertex_Index] != -1
                    && shortest_distance[End_Vertex_Index] != MAX_DISTANCE){
                        PathInformation* path = new PathInformation;
                        path->path_length = getPathWeightLength(first_shortest_path);
                        path->path_go_through_vertices = first_shortest_path;
                        top_k_shortest_path->insertLast(path);
                    }
                    // 此时说明找不到第一条最短路径 --> print outcome then delete all and return
                    else{
                        delete first_shortest_path;
                        printFinalOutcome(Top_K_Shortest_Path_Num,top_k_shortest_path);
                        delete[] visited;
                        delete[] shortest_distance;
                        delete[] pre_vertex_index;
                        delete[] transit_vertex_counter;
                        delete[] is_in_set;
                        delete top_k_shortest_path;
                        return;
                    }
                    // 既然已经找到了第一条最短路径, 最后我们找出剩余的k-1条路径即可
                        /*
                        实现难点 : 通过第一条最短路径生成候选路径 --> 根据候选路径的长度进行排序 --> 优先队列
                        */
                    VecList<PathInformation*>* candidatePaths = new VecList<PathInformation*>;
                    generateCandidatePaths(Start_Vertex,End_Vertex,top_k_shortest_path,candidatePaths);
                    int current_top_k_shortest_path_num = top_k_shortest_path->getLength();
                    while(current_top_k_shortest_path_num < Top_K_Shortest_Path_Num && !candidatePaths->isEmpty()){
                        deleteSamePath(top_k_shortest_path,candidatePaths);
                        PathInformation* new_shortest_path = candidatePaths->getEleAtPos(0);
                        candidatePaths->deleteEleAtPos(0);
                        if(new_shortest_path->path_go_through_vertices->getLength()-2<=max_transition_vertex_num){
                            top_k_shortest_path->insertLast(new_shortest_path);
                        }
                        generateCandidatePaths(Start_Vertex,End_Vertex,top_k_shortest_path,candidatePaths);
                        current_top_k_shortest_path_num++;
                    }
                    // 最后输出最终结果即可
                    printFinalOutcome(Top_K_Shortest_Path_Num,top_k_shortest_path);
                        cout << "The run time is: " <<(double)clock() / CLOCKS_PER_SEC << "s" << endl;
                    // 最后输出最终结果即可

                    // 最后清除内存
                    int actual_path_num = top_k_shortest_path->getLength();
                    for(int i=0;i<actual_path_num;i++){
                        PathInformation* now_path = top_k_shortest_path->getEleAtPos(i);
                        delete now_path->path_go_through_vertices;
                        delete now_path;
                    }
                    delete[] visited;
                    deleteMemory(shortest_distance, pre_vertex_index, is_in_set, transit_vertex_counter);
                    delete top_k_shortest_path;
                    for(int i=0;i<candidatePaths->getLength();i++){
                        PathInformation* now_path = candidatePaths->getEleAtPos(i);
                        delete now_path->path_go_through_vertices;
                        delete now_path;
                    }
                    delete candidatePaths;
                    // 最后清除内存
                }
            // ----------------------------------------  Paths ---------------------------------------- //

        //---------------------------------------------------------- 核心功能函数 ----------------------------------------------------------//

        //---------------------------------------------------------- 辅助函数 ----------------------------------------------------------//

            // ---------------------------------------- Others ---------------------------------------- //

                // ------------------------------ 初始化函数 ------------------------------//
                    bool* initVisitedArray(){
                        bool* visited = new bool[numVer];
                        for(int i=0;i<numVer;i++){
                            visited[i] = false;
                        }
                        return visited;
                    }

                    // int* initDistArray(){
                    //     int* dist = new int[numVer];
                    //     // 因为还不知道每个顶点的最短距离 --> set infinity
                    //     for(int i=0;i<numVer;i++){
                    //         dist[i] = MAX_DISTANCE;
                    //     }
                    //     return dist;
                    // }

                    // int* initPreArray(){
                    //     int* pre = new int[numVer];
                    //     // 因为还不知道每个顶点的上一个顶点的index --> set -1
                    //     for(int i=0;i<numVer;i++){
                    //         pre[i] = -1;
                    //     }
                    //     return pre;
                    // }

                    // int* initTransArray(){
                    //     int* trans = new int[numVer];
                    //     // 因为都还没有经过中继点 --> set -1
                    //     for(int i=0;i<numVer;i++){
                    //         trans[i] = -1;
                    //     }
                    //     return trans;
                    // }

                    // bool* initInSetArray(){
                    //     bool* isInS = new bool[numVer];
                    //     for(int i=0;i<numVer;i++){
                    //         isInS[i] = false;
                    //     }
                    //     return isInS;
                    // }

                    void initDijkstra(int* dist, int* pre, int* trans, bool* isInS){
                        for(int i=0;i<numVer;i++){
                            dist[i] = MAX_DISTANCE;
                            pre[i] = -1;
                            trans[i] = -1;
                            isInS[i] = false;
                        }
                    }
                // ------------------------------ 初始化函数 ------------------------------//

                // ------------------------------ 用于记录路径的结构体 ------------------------------//
                    struct PathInformation{
                        // 路径总权值长度
                        int path_length;
                        // 路径上经过的顶点
                        VecList<T>* path_go_through_vertices;
                    };
                // ------------------------------ 用于记录路径的结构体 ------------------------------//

                // ------------------------------ 记录删除边的结构体 ------------------------------//
                    struct RemoveEdge{
                        T Start_Vertex;
                        T End_Vertex;
                        int Delete_Edge_Weight;
                    };
                // ------------------------------ 记录删除边的结构体 ------------------------------//

            // ---------------------------------------- Others ---------------------------------------- //

            // ---------------------------------------- Ban Help ---------------------------------------- //

                // ------------------------------ 保存出边 ------------------------------//
                    void saveOutVertexAndEdgeWeights(int Ban_Vertex_Index, T Ban_Vertex){
                        for(int i=0;i<numVer;i++){
                            // 如果有出边
                            if(adjMatrix[Ban_Vertex_Index][i]!=0){
                                // 找出出顶点
                                T out_vertex = verList.getEleAtPos(i);
                                // 记录出边的权值
                                original_banned_vertex_edges[Ban_Vertex][out_vertex] = adjMatrix[Ban_Vertex_Index][i];
                                // 将出边禁用
                                adjMatrix[Ban_Vertex_Index][i] = 0;
                            }
                        }
                    }
                // ------------------------------ 保存出边 ------------------------------//

                // ------------------------------ 保存入边 ------------------------------//
                    void saveInVertexAndEdgeWeights(int Ban_Vertex_Index, T Ban_Vertex){
                        // cout << "Ban Vertex Index is " << Ban_Vertex_Index << endl;
                        for(int i=0;i<numVer;i++){
                            // 如果有入边
                            if(adjMatrix[i][Ban_Vertex_Index]!=0){
                                // 找出入顶点
                                T in_vertex = verList.getEleAtPos(i);
                                // 记录入边的权值
                                original_banned_vertex_edges[in_vertex][Ban_Vertex] = adjMatrix[i][Ban_Vertex_Index];
                                // 将入边禁用
                                adjMatrix[i][Ban_Vertex_Index] = 0;
                            }
                        }
                    }
                // ------------------------------ 保存入边 ------------------------------//

            // ---------------------------------------- Ban Help ---------------------------------------- //
            
            // ---------------------------------------- Unban Help ---------------------------------------- //

                // ------------------------------ Black List ------------------------------//
                    bool isVertexInBlackList(T Vertex){
                        // 实现难点 : find & end
                            // 如果 vertex在black list里面 <--> find() != end() <--> return true
                            // 如果 vertex在black list里面 <--> find() == end() <--> return false 
                        return black_list.find(Vertex) != black_list.end();
                    }
                // ------------------------------ Black List ------------------------------//

                // ------------------------------ 恢复顶点和边的权值 ------------------------------//
                    void recoverVertexAndEdgeWeight(T Unban_Vertex, int Unban_Vertex_Index){
                        /*
                        实现难点 : 关于unordered_map的使用
                            // it->first : 邻居节点
                            // it->second : 边的权值
                        */
                        for (typename unordered_map<T, int>::iterator it = original_banned_vertex_edges[Unban_Vertex].begin();it != original_banned_vertex_edges[Unban_Vertex].end(); ++it) {
                            // 获取邻居顶点和边权值
                            T vertex_neighbor = it->first;
                            int edge_weight = it->second;
                            // 找出解禁顶点的邻居的索引
                            int vertex_neighbor_index = verList.locateEle(vertex_neighbor);
                            // 如果邻居在Verlist里面 --> 恢复
                            if(vertex_neighbor_index != -1){
                                adjMatrix[Unban_Vertex_Index][vertex_neighbor_index] = edge_weight;
                                adjMatrix[vertex_neighbor_index][Unban_Vertex_Index] = edge_weight;
                            }
                        }
                    }
                // ------------------------------ 恢复顶点和边的权值 ------------------------------//
            // ---------------------------------------- Unban Help ---------------------------------------- //

            // ---------------------------------------- MaxTrans Help ---------------------------------------- //
            // ---------------------------------------- MaxTrans Help ---------------------------------------- //

            // ---------------------------------------- Path Help ---------------------------------------- //

                // ------------------------------ Dijkstra ------------------------------//
                    // ------------------------------ Dijkstra_init ------------------------------ //
                    void DijkstraInitDistAndPreNode(int sIndex, int* dist, int* preNode, bool* isInS, int* Transit_Vertex_Counter){
                        // traverse所有顶点
                        for(int i=0;i<numVer;i++){
                            // 如果顶点已经在集合里边了 --> continue
                            if(isInS[i]) continue;
                            // 此时保证了顶点不在集合里边 --> operate
                            /*
                            实现关键 : 跳过黑名单中的节点
                            */
                            T now_vertex = verList.getEleAtPos(i);
                            // 如果顶点在黑名单里边 --> 将该顶点的dis和prenode抹掉 并且 修改中继点个数 --> continue 
                            if(isVertexInBlackList(now_vertex)){
                                dist[i] = -1;
                                preNode[i] = -1;
                                Transit_Vertex_Counter[i] = -1;
                                continue;
                            }
                            // 此时说明顶点不在黑名单里边 --> modify该顶点的dis和prenode和transit counter
                            // 如果有直接连接的边
                                // 记录终止顶点到开始顶点的距离
                            dist[i] = (adjMatrix[sIndex][i]!=0 ? adjMatrix[sIndex][i] : -1);
                                // 记录前一个顶点的index
                            preNode[i] = (adjMatrix[sIndex][i]!=0 ? sIndex : -1);
                                // 记录中继点个数
                            Transit_Vertex_Counter[i] = (adjMatrix[sIndex][i]!=0 ? 0 : -1);
                        }
                    }
                    // ------------------------------ Dijkstra_init ------------------------------ //
                    void Dijkstra(T st, int* dist, int* preNode, bool* isInS, int* Transit_Vertex_Counter){
                        int sIndex = verList.locateEle(st);
                        dist[sIndex] = 0;
                        isInS[sIndex] = true;
                        // 起点是没有中继点的
                        Transit_Vertex_Counter[sIndex] = 0;
                        int countS = 1;
                        DijkstraInitDistAndPreNode(sIndex, dist, preNode, isInS, Transit_Vertex_Counter);
                            // printTransArr(Transit_Vertex_Counter);
                        while(countS < numVer){
                            int smallestDist = -1;
                            int corrIndex = -1;
                            for(int i=0;i<numVer;i++){
                                /*
                                实现关键 : 如果顶点在黑名单当中 --> 跳过顶点
                                */
                                T now_vertex = verList.getEleAtPos(i);
                                if(isVertexInBlackList(now_vertex)){
                                        // cout << "Skipping banned vertex: " << now_vertex << endl;
                                    continue;
                                } 
                                if(isInS[i]) continue;
                                if(dist[i]==-1) continue;
                                if(corrIndex==-1){
                                    smallestDist = dist[i];
                                    corrIndex = i;
                                }
                                else if(dist[i] < smallestDist){
                                    smallestDist = dist[i];
                                    corrIndex = i;
                                }
                            }
                            if(corrIndex==-1){
                                countS = numVer;
                                continue;
                            }
                            isInS[corrIndex] = true;
                            countS++;
                            for(int i=0;i<numVer;i++){
                                /*
                                实现关键 : 如果顶点在黑名单里面 --> 跳过该顶点
                                */
                                T now_vertex = verList.getEleAtPos(i);
                                if(isVertexInBlackList(now_vertex)) continue; 
                                if(isInS[i]) continue;
                                if(adjMatrix[corrIndex][i]==0) continue;
                                if(dist[corrIndex]==-1) continue;
                                int now_distance = dist[corrIndex] + adjMatrix[corrIndex][i];
                                // 起点sIndex 到 corrIndex 的中继点个数
                                int now_transition_vertex_counter = Transit_Vertex_Counter[corrIndex];
                                // 只要corrIndex不是起点 <--> 新增中继点 --> 那么now_transition_vertex_counter就等于Transit_Vertex_Counter[corrIndex]+1
                                if(corrIndex != sIndex){
                                    now_transition_vertex_counter +=1;
                                }
         
                                bool pass_transit_check = (max_transition_vertex_num<0) || (now_transition_vertex_counter<= max_transition_vertex_num);
                                bool can_update = (dist[i]== -1 || now_distance<dist[i]);
                                if(pass_transit_check && can_update){
                                    dist[i] = now_distance;
                                    preNode[i] = corrIndex;
                                    Transit_Vertex_Counter[i] = now_transition_vertex_counter;
                                }
                            }
                        }
                        // 此时已经完成了满足有限制的Dijkstra
                    }
                // ------------------------------ Dijkstra ------------------------------//


                // ------------------------------ 找出第一条最短路径 ------------------------------//
                    VecList<T>* findFirstShortestPath(T Start_Vertex,T End_Vertex,int* Shortest_Distance,int* Pre_Vertex_Index){
                        VecList<T>* path_go_through_vertices = new VecList<T>;
                        int end_vertex_index = verList.locateEle(End_Vertex);
                        // 如果终点不存在 --> return
                        if(end_vertex_index == -1){
                            return path_go_through_vertices;
                        }
                        // 如果终点不可达 --> return
                        if(Shortest_Distance[end_vertex_index] == -1){
                            return path_go_through_vertices;
                        }
                        // 如果终点不可达 --> return
                        if(Shortest_Distance[end_vertex_index] == MAX_DISTANCE){
                            return path_go_through_vertices;
                        }
                        // 此时说明终点是存在的 --> record path
                        int current_vertex_index = end_vertex_index;
                        // 不断地找出终点 until unreachable
                        while(current_vertex_index !=-1){
                            T current_vertex = verList.getEleAtPos(current_vertex_index);
                            path_go_through_vertices->insertLast(current_vertex);
                            current_vertex_index = Pre_Vertex_Index[current_vertex_index];
                        }
                        // 此时已经放完了(但是因为我们从后往前放的) --> 我们还要翻转一下这个向量表 --> return VList
                        path_go_through_vertices->reverseEleInVecList();
                        // 此时已经得到结果 --> return result
                        return path_go_through_vertices;
                    }
                // ------------------------------ 找出第一条最短路径 ------------------------------//

                // ------------------------------ 找出顶点集的路径权值和 ------------------------------//
                    int getPathWeightLength(VecList<T>* Path_Go_Through_Vertices){
                        // 如果只有一个顶点 --> return 0
                        if(Path_Go_Through_Vertices->getLength() == 0) return 0;
                        int weight_length = 0;
                        // 遍历路径中的每一对相邻顶点，计算它们之间的边权值
                        for (int i = 0; i < Path_Go_Through_Vertices->getLength() - 1; i++) {
                            T start_vertex = Path_Go_Through_Vertices->getEleAtPos(i);
                            T end_vertex = Path_Go_Through_Vertices->getEleAtPos(i + 1);
                            
                            // 获取当前顶点和下一个顶点的索引
                            int start_vertex_index = verList.locateEle(start_vertex);
                            int end_vertex_index = verList.locateEle(end_vertex);
                            
                            // 根据邻接矩阵找到两个顶点之间的边权值
                            int edge_weight = adjMatrix[start_vertex_index][end_vertex_index];
                            
                            // 累加边的权值
                            weight_length += edge_weight;
                        }
                        
                        return weight_length;
                    }
                // ------------------------------ 找出顶点集的路径权值和 ------------------------------//

                // ------------------------------ 基于候选路径的权值长度插入到候选路径向量表 ------------------------------//
                    void insertCandidatePathBaseOnWeight(VecList<PathInformation*>* Candidates_Paths, PathInformation* New_Candidate_Path){
                        int now_path_num = Candidates_Paths->getLength();
                        for(int i=0;i<now_path_num;i++){
                            int new_candidate_path_length = New_Candidate_Path->path_length;
                            PathInformation* origin_candidates_path = Candidates_Paths->getEleAtPos(i);
                            int origin_candidates_path_length = origin_candidates_path->path_length;
                            // 只要新的路径比原来的路径长度要短 --> insert and return
                            if(new_candidate_path_length < origin_candidates_path_length){
                                Candidates_Paths->insertEleAtPos(i,New_Candidate_Path);
                                return;
                            }
                        }
                        // 此时新的路径比原来的路径长度要长 --> 插入到队尾
                        Candidates_Paths->insertLast(New_Candidate_Path);
                    }
                // ------------------------------ 基于候选路径的权值长度插入到候选路径向量表 ------------------------------//

                // ------------------------------ 生成候选路径 ------------------------------//
                    void generateCandidatePaths(T Start_Vertex, T End_Vertex, VecList<PathInformation*>* Top_K_Shortest_Paths, VecList<PathInformation*>* Candidates_Paths){
                        int end_vertex_index = verList.locateEle(End_Vertex);
                        int now_top_k_shortest_path_num = Top_K_Shortest_Paths->getLength();
                        PathInformation* last_shortest_path = Top_K_Shortest_Paths->getEleAtPos(now_top_k_shortest_path_num-1);
                        VecList<T>* last_shortest_path_go_through_vertex = last_shortest_path->path_go_through_vertices;
                        int last_shortest_path_go_through_vertex_num = last_shortest_path_go_through_vertex->getLength();
                        VecList<RemoveEdge*>* remove_edges_set = new VecList<RemoveEdge*>;
                        addEdgesUsedToBeDeleted(remove_edges_set,Top_K_Shortest_Paths);
                        // 从头开始, 以每个节点为偏转节点, 不断生成新的候选路径
                        for(int i=0;i<last_shortest_path_go_through_vertex_num-1;i++){
                            VecList<T>* root_path = new VecList<T>;
                            addVertexToRootPath(root_path, last_shortest_path_go_through_vertex,i+1);
                            // 获取偏转节点
                            T deviation_node = last_shortest_path_go_through_vertex->getEleAtPos(i);
                            int deviation_node_index = verList.locateEle(deviation_node);
                            T deviation_next_node = last_shortest_path_go_through_vertex->getEleAtPos(i+1);
                            int deviation_next_node_index = verList.locateEle(deviation_next_node);
                            RemoveEdge* new_remove_edge = new RemoveEdge;
                            new_remove_edge->Start_Vertex = deviation_node;
                            new_remove_edge->End_Vertex = deviation_next_node;
                            new_remove_edge->Delete_Edge_Weight = adjMatrix[deviation_node_index][deviation_next_node_index];
                            // 如果找不到 --> 插入
                            removeEdgesWithSamePrefix(root_path, Top_K_Shortest_Paths);
                            // 重新执行Dijkstra
                            int* new_shortest_distance = new int[numVer];
                            int* new_pre_vertex_index = new int[numVer];    
                            bool* new_is_in_set = new bool[numVer];
                            int* new_transit_vertex_counter = new int[numVer];
                            initDijkstra(new_shortest_distance,new_pre_vertex_index,new_transit_vertex_counter,new_is_in_set)
                            Dijkstra(deviation_node, new_shortest_distance, new_pre_vertex_index, new_is_in_set, new_transit_vertex_counter);
                            // 如果起点去不到终点 --> delete memory and continue
                            if(new_shortest_distance[end_vertex_index] == -1){
                                delete root_path;
                                deleteMemory(new_shortest_distance, new_pre_vertex_index, new_is_in_set, new_transit_vertex_counter);
                                continue;
                            }
                            // 此时起点可以去到终点
                            VecList<T>* generate_new_shortest_path = findFirstShortestPath(deviation_node, End_Vertex, new_shortest_distance, new_pre_vertex_index);
                            // 如果有重复出现的顶点 --> 跳过
                            // 此时说明没有重复出现的顶点 --> 继续搜索
                            // 同样地, 如果找的到新的最短路径 --> record
                                /*
                                实现关键 : root path中的节点不变 , new shortrest path中的节点改变
                                */
                            int root_path_length = root_path->getLength();
                            int generate_new_shortest_path_length = generate_new_shortest_path->getLength();
                            if(root_path_length > 0){
                                VecList<T>* new_shortest_path = new VecList<T>;
                                for(int i=0;i<root_path_length;i++){
                                    T root_path_vertex = root_path->getEleAtPos(i);
                                    new_shortest_path->insertLast(root_path_vertex);
                                }
                                /*
                                实现关键 : 因为index = 0 已经入队 --> 只需从index = 1开始即可
                                */
                                for(int i=1;i<generate_new_shortest_path_length;i++){
                                    T generate_new_shortest_path_vertex = generate_new_shortest_path->getEleAtPos(i);
                                    new_shortest_path->insertLast(generate_new_shortest_path_vertex);
                                }  
                                // recoverOriginalGraph(original_graph);
                                // 此时将新生成的路径放入到候选路径里边
                                PathInformation* new_candidates_path = new PathInformation;
                                new_candidates_path->path_length = getPathWeightLength(generate_new_shortest_path) + addPathWeightDeleted(root_path, remove_edges_set);
                                new_candidates_path->path_go_through_vertices = new_shortest_path;
                                insertCandidatePathBaseOnWeight(Candidates_Paths, new_candidates_path);
                            }
                            // 此时已经找完了一条路径生成的所有候选路径 --> delete memory
                            delete root_path;
                            delete generate_new_shortest_path;
                            deleteMemory(new_shortest_distance, new_pre_vertex_index, new_is_in_set, new_transit_vertex_counter);
                        }
                        // 此时已经通过一条路径生成了新的候选路径
                        recoverOriginalGraph(remove_edges_set);
                    }
                // ------------------------------ 生成候选路径 ------------------------------//

            // ---------------------------------------- Path Help ---------------------------------------- //

                // ------------------------------ 移除具有共享前缀的边 ------------------------------//
                    void removeEdgesWithSamePrefix(VecList<T>* root_path, VecList<PathInformation*>* Top_K_Shortest_Paths) {
                        for (int i = 0; i < Top_K_Shortest_Paths->getLength(); i++) {
                            PathInformation* tmp_path = Top_K_Shortest_Paths->getEleAtPos(i);
                            VecList<T>* tmp_path_go_through_vertices = tmp_path->path_go_through_vertices;
                            bool has_same_prefix = true;
                            // 检查前缀是否相同
                            for (int j = 0; j < root_path->getLength(); j++) {
                                T root_path_vertex = root_path->getEleAtPos(j) ;
                                T tmp_path_vertex = tmp_path_go_through_vertices->getEleAtPos(j);
                                // 如果没有共同的前缀 --> set false and break
                                if (root_path_vertex != tmp_path_vertex) {
                                    has_same_prefix = false;
                                    break;
                                }
                            }
                            // 此时有共同的前缀 --> 移除前缀后的第一条边
                            if (has_same_prefix && tmp_path_go_through_vertices->getLength() > root_path->getLength()) {
                                T start_vertex = tmp_path_go_through_vertices->getEleAtPos(root_path->getLength() - 1);
                                T end_vertex = tmp_path_go_through_vertices->getEleAtPos(root_path->getLength());
                                int start_index = verList.locateEle(start_vertex);
                                int end_index = verList.locateEle(end_vertex);
                                adjMatrix[start_index][end_index] = 0;
                                adjMatrix[end_index][start_index] = 0;
                            }
                        }
                    }
                // ------------------------------ 移除具有共享前缀的边 ------------------------------//

                // ------------------------------ 基于删除的边修改邻接矩阵 ------------------------------//
                    void ModifyGraphBaseOnRemoveEdges(VecList<RemoveEdge*>* Remove_Edges_Set){
                        // 如果是空的 --> return
                        if(Remove_Edges_Set->isEmpty()){
                            return;
                        }
                        // 此时不是空的 --> modify
                        for(int i=0;i<Remove_Edges_Set->getLength();i++){
                            RemoveEdge* tmp_remove_edge = Remove_Edges_Set->getEleAtPos(i);
                            T tmp_remove_edge_start_vertex = tmp_remove_edge->Start_Vertex;
                            int tmp_remove_edge_start_vertex_index = verList.locateEle(tmp_remove_edge_start_vertex);
                            T tmp_remove_edge_end_vertex = tmp_remove_edge->End_Vertex;
                            int tmp_remove_edge_end_vertex_index = verList.locateEle(tmp_remove_edge_end_vertex);
                            adjMatrix[tmp_remove_edge_start_vertex_index][tmp_remove_edge_end_vertex_index] = 0;
                            adjMatrix[tmp_remove_edge_end_vertex_index][tmp_remove_edge_start_vertex_index] = 0;                 
                        }
                    }
                // ------------------------------ 基于删除的边修改邻接矩阵 ------------------------------//
                // ------------------------------ 添加曾经删除过的边 ------------------------------//
                    void addEdgesUsedToBeDeleted(VecList<RemoveEdge*>* Remove_Edges_Set, VecList<PathInformation*>* Top_K_Shortest_Paths){
                        // 如果是空 --> return
                        if(Top_K_Shortest_Paths->isEmpty()){
                            return;
                        }
                        // 此时不是空 --> 可以删除
                        for(int i=0;i<Top_K_Shortest_Paths->getLength();i++){
                            PathInformation* tmp_path_info = Top_K_Shortest_Paths->getEleAtPos(i);
                            VecList<T>* tmp_path_go_through_vertex = tmp_path_info->path_go_through_vertices;
                            for(int j=0;j<tmp_path_go_through_vertex->getLength()-1;j++){
                                RemoveEdge* remove_edge = new RemoveEdge;
                                T start_vertex = tmp_path_go_through_vertex->getEleAtPos(j);
                                T end_vertex = tmp_path_go_through_vertex->getEleAtPos(j+1);
                                int start_vertex_index = verList.locateEle(start_vertex);
                                int end_vertex_index = verList.locateEle(end_vertex);
                                remove_edge->Start_Vertex = start_vertex;
                                remove_edge->End_Vertex = end_vertex;
                                remove_edge->Delete_Edge_Weight = adjMatrix[start_vertex_index][end_vertex_index];
                                // 如果已经存在了 --> delete and continue
                                if(isEdgeAlreadyRemoved(remove_edge,Remove_Edges_Set)){
                                    delete remove_edge;
                                    continue;
                                }
                                // 如果还没存在 --> insert last
                                else{
                                    Remove_Edges_Set->insertLast(remove_edge);
                                }
                            }
                        }
                    }
                // ------------------------------ 添加曾经删除过的边 ------------------------------//
                // ------------------------------ 找出已经删除的边的所有权值 ------------------------------//
                    int addPathWeightDeleted(VecList<T>* Root_Path, VecList<RemoveEdge*>* Remove_Edge_Set){
                        // 如果还没删除边 --> 删除边的权值为0
                        if(Remove_Edge_Set->isEmpty()){
                            return 0;
                        }
                        // 此时有删除边
                        else{
                            int delete_edge_weight = 0;
                            for(int i=0;i<Root_Path->getLength()-1;i++){
                                T root_path_start_vertex = Root_Path->getEleAtPos(i);
                                T root_path_end_vertex = Root_Path->getEleAtPos(i+1);
                                for(int j=0;j<Remove_Edge_Set->getLength();j++){
                                    RemoveEdge* tmp_remove_edge = Remove_Edge_Set->getEleAtPos(j);
                                    T start_vertex = tmp_remove_edge->Start_Vertex;
                                    T end_vertex = tmp_remove_edge->End_Vertex;
                                    if((start_vertex == root_path_start_vertex && end_vertex == root_path_end_vertex)
                                    ||(end_vertex == root_path_start_vertex && start_vertex == root_path_end_vertex)){
                                        delete_edge_weight += tmp_remove_edge->Delete_Edge_Weight;
                                    }
                                }
                            }
                            return delete_edge_weight;
                        }
                    }
                // ------------------------------ 找出已经删除的边的所有权值 ------------------------------//

                // ------------------------------ 可以找到已经删除过的边 ------------------------------//
                    bool isEdgeAlreadyRemoved(RemoveEdge* Remove_Edge, VecList<RemoveEdge*>* Remove_Edge_Set){
                        for(int i=0;i<Remove_Edge_Set->getLength();i++){
                            RemoveEdge* tmp_remove_edge = Remove_Edge_Set->getEleAtPos(i);
                            if(((Remove_Edge->Start_Vertex == tmp_remove_edge->Start_Vertex && Remove_Edge->End_Vertex == tmp_remove_edge->End_Vertex)
                            || (Remove_Edge->Start_Vertex == tmp_remove_edge->End_Vertex && Remove_Edge->End_Vertex == tmp_remove_edge->Start_Vertex)) 
                            && (Remove_Edge->Delete_Edge_Weight == tmp_remove_edge->Delete_Edge_Weight)){
                                return true;
                            }
                        }
                        return false;
                    }
                // ------------------------------ 可以找到已经删除过的边 ------------------------------//
                // ------------------------------ 判断生成的路径是否已经有之前的顶点 ------------------------------//
                    bool hasSomeVertexhasBeenVisited(VecList<T>* Root_Path, VecList<T>* Generate_New_Path){
                        for(int i=1;i<Generate_New_Path->getLength();i++){
                            T Generate_New_Path_vertex = Generate_New_Path->getEleAtPos(i);
                            // 如果找到了之前已经出现过的顶点 --> return true
                            // 如果找不到 --> 一直往后边找
                            if(Root_Path->locateEle(Generate_New_Path_vertex) == -1){
                                continue;
                            }
                            // 如果找得到 --> 说明有重复顶点 --> return true
                            if(Root_Path->locateEle(Generate_New_Path_vertex) >= 0){
                                return true;
                            }
                        }
                        // 此时说明新的路径中没有之前出现过的顶点 --> return false
                        return false;
                    }
                // ------------------------------ 判断生成的路径是否已经有之前的顶点 ------------------------------//

                // ------------------------------ 移除已经走过的路 ------------------------------//
                    void removePathsHaveBeenThrough(VecList<T>* Root_Path){
                        if(Root_Path->getLength() <= 1) return;
                        for(int i=0;i<Root_Path->getLength()-1;i++){
                            T start_vertex = Root_Path->getEleAtPos(i);
                            T end_vertex = Root_Path->getEleAtPos(i+1);
                            int start_vertex_index = verList.locateEle(start_vertex);
                            int end_vertex_index = verList.locateEle(end_vertex);
                            adjMatrix[start_vertex_index][end_vertex_index] = 0;
                            adjMatrix[end_vertex_index][start_vertex_index] = 0;
                        }
                    }
                // ------------------------------ 移除已经走过的路 ------------------------------//

                // ------------------------------ 判断是否为相同路径 ------------------------------//
                    bool isSamePaths(VecList<T>* Path_1, VecList<T>* Path_2){
                        // 如果顶点个数不同 --> 肯定不是
                        if(Path_1->getLength() != Path_2->getLength()) return false;
                        // 此时顶点个数相同 --> 判断是不是经过路径上每个顶点都一样
                            // 如果是 --> 相同路径
                            // 如果不是 --> 不是相同路径
                        for(int i=0;i<Path_1->getLength();i++){
                            T path1_vertex = Path_1->getEleAtPos(i);
                            T path2_vertex = Path_2->getEleAtPos(i);
                            if(path1_vertex != path2_vertex){
                                return false;
                            }
                        }
                        return true;
                    }
                // ------------------------------ 判断是否为相同路径 ------------------------------//

                // ------------------------------ 删除Candidates path中的重复路径(这些重复路径在TopK也出现过) ------------------------------//
                    void deleteSamePath(VecList<PathInformation*>* Top_K_Shortest_Paths, VecList<PathInformation*>* Candidates_Paths){
                        // 如果为空 --> 不做任何操作
                        if(Top_K_Shortest_Paths->isEmpty()) return;
                        if(Candidates_Paths->isEmpty()) return;
                        // 此时保证了不为空 --> 删除候选路径中的重复项
                        for(int i=0;i<Top_K_Shortest_Paths->getLength();i++){
                            PathInformation* tmp_path = Top_K_Shortest_Paths->getEleAtPos(i);
                            VecList<T>* tmp_path_go_through_vertex = tmp_path->path_go_through_vertices;
                            for(int j =0;j<Candidates_Paths->getLength();){
                                PathInformation* now_cadidate_path = Candidates_Paths->getEleAtPos(j);
                                VecList<T>* now_cadidate_path_go_through_vertex = now_cadidate_path->path_go_through_vertices;
                                // 如果是相同路径 --> 进行删除
                                if(isSamePaths(tmp_path_go_through_vertex, now_cadidate_path_go_through_vertex)){
                                    delete now_cadidate_path->path_go_through_vertices;
                                    delete now_cadidate_path;
                                    Candidates_Paths->deleteEleAtPos(j);
                                    /*
                                    教训 : 之前也是这么错的
                                    */
                                }
                                // 如果不是相同路径 --> 继续往后进行比较
                                else{
                                    j++;
                                }
                            }
                        }
                    }
                // ------------------------------ 删除Candidates path中的重复路径(这些重复路径在TopK也出现过) ------------------------------//

                // ------------------------------ 移除root path中的路径 ------------------------------//
                    void removeEdgesBaseOnRootPath(int Root_Path_Length, VecList<T>* Last_Shortest_Path_Go_Through_Vertex){
                        for(int i=Root_Path_Length-1;i<Last_Shortest_Path_Go_Through_Vertex->getLength()-1;i++){
                            T start_vertex_in_root_path = Last_Shortest_Path_Go_Through_Vertex->getEleAtPos(i);
                            T end_vertex_in_root_path = Last_Shortest_Path_Go_Through_Vertex->getEleAtPos(i+1);
                            int start_vertex_in_root_path_index = verList.locateEle(start_vertex_in_root_path);
                            int end_vertex_in_root_path_index = verList.locateEle(end_vertex_in_root_path);
                            adjMatrix[start_vertex_in_root_path_index][end_vertex_in_root_path_index] = 0;
                            adjMatrix[end_vertex_in_root_path_index][start_vertex_in_root_path_index] = 0;
                        }
                    }
                // ------------------------------ 移除root path中的路径 ------------------------------//

                // ------------------------------ 释放内存 ------------------------------//
                    void deleteMemory(int* shortest_distance,int* pre_vertex_index,bool* is_in_set,int* transit_vertex_counter){
                        delete[] shortest_distance;
                        delete[] pre_vertex_index;
                        delete[] is_in_set;
                        delete[] transit_vertex_counter;
                    }
                // ------------------------------ 释放内存 ------------------------------//


                // ------------------------------ 记录原始图 ------------------------------//
                    int** backUpOriginalGraph(){
                        // 初始化备份原始图
                        int** back_up_graph = new int*[numVer];
                        for(int i=0;i<numVer;i++){
                            back_up_graph[i] = new int[numVer];
                        }
                        // 记录原始图
                        for(int i=0;i<numVer;i++){
                            for(int j=0;j<numVer;j++){
                                back_up_graph[i][j] = adjMatrix[i][j];
                            }
                        }
                        return back_up_graph;
                    }
                // ------------------------------ 记录原始图 ------------------------------//

                // ------------------------------ 恢复原始图 ------------------------------//
                    void recoverOriginalGraph(VecList<RemoveEdge*>* Remove_Edges_Set){
                        // 如果是空的 --> return
                        if(Remove_Edges_Set->isEmpty()){
                            return;
                        }
                        // 如果不是空的 --> recover
                        else{
                            for(int i=0;i<Remove_Edges_Set->getLength();i++){
                                RemoveEdge* tmp_remove_edge = Remove_Edges_Set->getEleAtPos(i);
                                T start_vertex = tmp_remove_edge->Start_Vertex;
                                T end_vertex = tmp_remove_edge->End_Vertex;
                                int start_vertex_index = verList.locateEle(start_vertex);
                                int end_vertex_index = verList.locateEle(end_vertex);
                                int edge_weight = tmp_remove_edge->Delete_Edge_Weight;
                                adjMatrix[start_vertex_index][end_vertex_index] = edge_weight;
                                adjMatrix[end_vertex_index][start_vertex_index] = edge_weight;
                            }
                        }
                        
                    }
                // ------------------------------ 恢复原始图 ------------------------------//

                // ------------------------------ 将原始顶点加入到root path ------------------------------//
                    void addVertexToRootPath(VecList<T>* new_vclist,VecList<T>* old_vclist,int length){
                        // printVecList(new_vclist);
                        for(int i=0;i<length;i++){
                            T to_set_val = old_vclist->getEleAtPos(i);
                            new_vclist->insertLast(to_set_val);
                        }
                    }
                // ------------------------------ 将原始顶点加入到root path ------------------------------//

                // ------------------------------ Outcome ------------------------------//
                    void printFinalOutcome(int Top_K_Shortest_Path_Num, VecList<PathInformation*>* Top_K_Shortest_Path){
                        // printVertexInBlackList();
                        // printMaxTransitionVertices();
                        cout << endl;
                        printTopKShortestPath(Top_K_Shortest_Path_Num,Top_K_Shortest_Path);
                    }

                    void printTopKShortestPath(int Top_K_Shortest_Path_Num, VecList<PathInformation*>* Top_K_Shortest_Path){
                        cout << "Top " << Top_K_Shortest_Path_Num << " Shortest Paths." << endl;
                        int actual_shortest_path_num = Top_K_Shortest_Path->getLength();
                        cout << "There are " << actual_shortest_path_num << " paths." << endl;;
                        for(int i=0;i < actual_shortest_path_num;i++){
                            PathInformation* shortest_path = Top_K_Shortest_Path->getEleAtPos(i);
                            VecList<T>* shortest_path_go_through_vertices = shortest_path->path_go_through_vertices;
                                // printVecList(shortest_path_go_through_vertices);
                            cout << "Path " << i << ": ";
                            int shortest_path_go_through_vertices_num = shortest_path_go_through_vertices->getLength();
                                // cout << "Shortest path go through vertices num is " << shortest_path_go_through_vertices_num << endl;
                            for(int j=0;j<shortest_path_go_through_vertices_num;j++){
                                cout << shortest_path_go_through_vertices->getEleAtPos(j) ;
                                    // cout << "j is " << j << endl;
                                    // cout << (j==shortest_path_go_through_vertices_num ? " j = shortest_path_go_through_vertices_num": "j != shortest_path_go_through_vertices_num") << endl;
                                cout << (j!=shortest_path_go_through_vertices_num-1 ?  " -> " : " ");
                            }
                            cout << ", distance = " << shortest_path->path_length << endl;
                        }
                    }
                // ------------------------------ Outcome ------------------------------//
            // ---------------------------------------- Path Help ---------------------------------------- //

        //---------------------------------------------------------- 辅助函数 ----------------------------------------------------------//

        //---------------------------------------------------------- 测试函数 ----------------------------------------------------------//

            // ------------------------------ Black List ------------------------------//
                void printVertexInBlackList(){
                    cout << "Vertices in Black List are: ";
                    for(T vertex : black_list){
                        cout << vertex << " ";
                    }
                    cout << endl;
                }
            // ------------------------------ Black List ------------------------------//

            // ------------------------------ Max Transition Vertices ------------------------------//
                void printMaxTransitionVertices(){
                    cout << "Max transition vertices number is set: " << max_transition_vertex_num << endl;
                }
            // ------------------------------ Max Transition Vertices ------------------------------//

            // ------------------------------ 打印邻接矩阵 ------------------------------//
                void printAdjencyMatrix(){
                    for(int i=0;i<numVer;i++){
                        for(int j=0;j<numVer;j++){
                            cout << adjMatrix[i][j] << " ";
                        }
                        cout << endl;
                    }
                    cout << endl;
                }
            // ------------------------------ 打印邻接矩阵 ------------------------------//

            // ------------------------------ 打印顶点集 ------------------------------//
                void printVerList(){
                    cout << "Verlist: " ;
                    for(int i=0;i<numVer;i++){
                        cout << verList.getEleAtPos(i) << " ";
                    }
                    cout << endl;
                }
            // ------------------------------ 打印顶点集 ------------------------------//

            // ------------------------------ 打印向量表 ------------------------------//
                void printVecList(VecList<T>* vclist){
                    if(vclist->getLength()==0){
                        cout << "Veclist is NULL!";
                        return;
                    }
                    cout << "Veclist: ";
                    for(int i=0;i<vclist->getLength();i++){
                        cout << vclist->getEleAtPos(i) << " ";
                    }
                    cout << endl;
                }
            // ------------------------------ 打印向量表 ------------------------------//

            // ------------------------------ 打印最短距离集 ------------------------------//
                void printDistArray(T Start_Vertex,int* dist){
                    for(int i=0;i<numVer;i++){
                        cout << "Frome Vertex " << Start_Vertex << " to Vertex " << verList.getEleAtPos(i) << "'s shortest distance is " << dist[i] << endl;
                    }
                    cout << endl;
                }
            // ------------------------------ 打印最短距离集 ------------------------------//

            // ------------------------------ 打印maxTrans集 ------------------------------//
                void printTransArr(int* Transit_Vertex_Counter){
                    for(int i=0;i<numVer;i++){
                            cout << "Vertex " << verList.getEleAtPos(i) << " 's transit num is " << Transit_Vertex_Counter[i] << endl;
                    }
                    cout << endl;
                }
            // ------------------------------ 打印maxTrans集 ------------------------------//

            // ------------------------------ 打印前K条最短路径 ------------------------------//
                void testTopKShortestPath(int Top_K_Shortest_Paths_Num, VecList<PathInformation*>* Top_K_Shortest_Path){
                    int actual_shortest_path_num = Top_K_Shortest_Path->getLength();
                    cout << "There are " << actual_shortest_path_num << " paths." << endl;;
                    for(int i=0;i < actual_shortest_path_num;i++){
                        PathInformation* shortest_path = Top_K_Shortest_Path->getEleAtPos(i);
                        VecList<T>* shortest_path_go_through_vertices = shortest_path->path_go_through_vertices;
                        cout << "Path " << i << ": ";
                        int shortest_path_go_through_vertices_num;
                        for(int j=0;j<shortest_path_go_through_vertices_num;j++){
                            cout << shortest_path_go_through_vertices->getEleAtPos(j) ;
                            cout << (j!=shortest_path_go_through_vertices_num ?  " -> " : " ");
                        }
                        cout << ", distance = " << shortest_path->path_length << endl;
                    }
                }
            // ------------------------------ 打印前K条最短路径 ------------------------------//

            // ------------------------------ 打印路径信息 ------------------------------//
                void printPathInformation(PathInformation* Path_Info){
                    cout << "Path length is " << Path_Info->path_length << endl;
                    printVecList(Path_Info->path_go_through_vertices);
                }
            // ------------------------------ 打印路径信息 ------------------------------//

            // ------------------------------ 打印候选路径 ------------------------------//
                void printCandidatesPath(VecList<PathInformation*>* Candidates_Paths){
                    for(int i=0;i<Candidates_Paths->getLength();i++){
                        printPathInformation(Candidates_Paths->getEleAtPos(i));
                    }
                    cout << endl;
                }
            // ------------------------------ 打印路径信息 ------------------------------//

            // ------------------------------ 打印移除路径向量表信息 ------------------------------//
                void printRemoveEdgesVecInfo(VecList<RemoveEdge*>* Remove_Edges_Set){
                    // 如果为空 --> 啥都没有
                    if(Remove_Edges_Set->isEmpty()){
                        cout << "Remove Edges Set is Empty!" << endl;
                    }
                    // 如果不为空 --> cout
                    else{
                        for(int i=0;i<Remove_Edges_Set->getLength();i++){
                            RemoveEdge* tmp_remove_edge = Remove_Edges_Set->getEleAtPos(i);
                            T tmp_start_vertex = tmp_remove_edge->Start_Vertex;
                            T tmp_end_vertex = tmp_remove_edge->End_Vertex;
                            int tmp_remove_edge_weight = tmp_remove_edge->Delete_Edge_Weight;
                            cout << "Start Vertex is " << tmp_start_vertex << endl;
                            cout << "End Vertex is " << tmp_end_vertex << endl;
                            cout << "Remove Edge Weight is " << tmp_remove_edge_weight << endl;
                            cout << endl;
                        }
                    }
                }
            // ------------------------------ 打印移除路径向量表信息 ------------------------------//

            // ------------------------------ 打印移除路径信息 ------------------------------//
                void printRemoveEdgesInfo(RemoveEdge* Remove_Edges){
                    T tmp_start_vertex = Remove_Edges->Start_Vertex;
                    T tmp_end_vertex = Remove_Edges->End_Vertex;
                    int tmp_remove_edge_weight = Remove_Edges->Delete_Edge_Weight;
                    cout << "Start Vertex is " << tmp_start_vertex << endl;
                    cout << "End Vertex is " << tmp_end_vertex << endl;
                    cout << "Remove Edge Weight is " << tmp_remove_edge_weight << endl;
                    cout << endl;
                }
            // ------------------------------ 打印移除路径信息 ------------------------------//

        //---------------------------------------------------------- 测试函数 ----------------------------------------------------------//
};

// ---------------------------------------------------------- 初始化函数 ---------------------------------------------------------- //
    template<class T>
    T* addVertices(int vertex_num){
        T* vertex_list = new T[vertex_num];
        for(int i=0;i<vertex_num;i++){
            T vertex_name;
            cin >> vertex_name;
            vertex_list[i] = vertex_name;
        }
        return vertex_list;
    }

    template<class T>
    void addEdge(AMGraph<T>* new_graph,int edge_num){
        for(int i=0;i<edge_num;i++){
            T start_vertex;
            T end_vertex;
            int edge_weight;
            cin >> start_vertex >> end_vertex >> edge_weight;
            new_graph->addEdge(start_vertex,end_vertex,edge_weight);
        }
    }

    template<class T>
    AMGraph<T>* initNewGraph(int vertex_num,int edge_num){
        T* vertex_list = addVertices<T>(vertex_num);
        AMGraph<T>* new_graph = new AMGraph<T>(vertex_list,vertex_num);
        // 此时已经用完了vertex list --> delete
        delete vertex_list;
        addEdge(new_graph,edge_num);
        return new_graph;
    }
// ---------------------------------------------------------- 初始化函数 ---------------------------------------------------------- //

// ---------------------------------------------------------- Solve ---------------------------------------------------------- //
    template<class T>
    void Solve(string command,AMGraph<T>* Graph){
        if(command == "ban"){
            T ban_vertex;
            cin >> ban_vertex;
            Graph->addVertexIntoBlackList(ban_vertex);
            Graph->printVertexInBlackList();
        }
        else if(command == "unban"){
            T unban_vertex;
            cin >> unban_vertex;
            Graph->removeVertexOutOfBlackList(unban_vertex);
            Graph->printVertexInBlackList();
        }
        else if(command == "maxTrans"){
            int max_trans_num;
            cin >> max_trans_num;
            Graph->limitMaxTransitVertexAs(max_trans_num);
        }
        else if(command == "paths"){
            T start_vertex;
            T end_vertex;
            int top_k_shortest_path;
            cin >> start_vertex >> end_vertex >> top_k_shortest_path;
            Graph->findTopKShortestPath(start_vertex,end_vertex,top_k_shortest_path);
            // Graph->printFinalOutcome();
        }
        else{
            cout <<  "Invalid command!" << endl;
        }
    }
// ---------------------------------------------------------- Solve ---------------------------------------------------------- //

int main(){
        // test();
    // ------------------------------ 初始化 ------------------------------//
    int m,n;
    cin >> m >> n;
    AMGraph<char>* Graph = initNewGraph<char>(m,n);
    // Graph->printGraph();
    // ------------------------------ 初始化 ------------------------------//
    // ------------------------------ 交互页面 ------------------------------//
    string command;
    while(true){
        cout << endl;
        cout << "Please enter your next request." << endl;
        cout << "Options include: ban X, unban X, maxTrans Y, paths A B k, quit." << endl;
        cout << endl;
        cin >> command;
        if(command == "quit"){
            break;
        }
        else{
            Solve(command,Graph);
        }
    }
    cout << "Program ended!" << endl;
    // ------------------------------ 交互页面 ------------------------------//
}