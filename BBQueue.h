#pragma once
#include "structures_2.h"
#include "BBNode.h"
#include <queue>
#include <stack>

using namespace std;

class BBQueue
{
public:

	BBQueue();
	~BBQueue();

	void setStrategy(BranchingStrategy strategy);

	void push(BBNode element);
	
	BBNode pop();

	long long int size();

	bool empty();

	long long int processed();

	void incrementProcessed();

private:

	queue<BBNode> bfsQueue_;
	priority_queue<BBNode, vector<BBNode>, BBNodeComparer> bfsPriorityQueue_;
	stack<BBNode> dfsStack_;

	BranchingStrategy strategy_;

	long long int processed_;
};

