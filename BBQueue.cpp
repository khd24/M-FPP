#include "BBQueue.h"

BBQueue::BBQueue() {
	strategy_ = BranchingStrategy::BstFS;
	processed_ = 0;
}

BBQueue::~BBQueue() {

}

void BBQueue::setStrategy(BranchingStrategy strategy) {
	strategy_ = strategy;
}

void BBQueue::push(BBNode element) {
	switch (strategy_) {
	case BranchingStrategy::BFS:
		return bfsQueue_.push(element);
	case BranchingStrategy::BstFS:
		return bfsPriorityQueue_.push(element);
	case BranchingStrategy::DFS:
		return dfsStack_.push(element);
	}
}

BBNode BBQueue::pop() {
	
	BBNode ret;

	switch (strategy_) {
	case BranchingStrategy::BFS:
		ret = bfsQueue_.front();
		bfsQueue_.pop();
		break;
	case BranchingStrategy::BstFS:
		ret = bfsPriorityQueue_.top();
		bfsPriorityQueue_.pop();
		break;
	case BranchingStrategy::DFS:
		ret = dfsStack_.top();
		dfsStack_.pop();
		break;
	default:
		ret = bfsPriorityQueue_.top();
		bfsPriorityQueue_.pop();
		break;
	}

	++processed_;

	return ret;
}

long long int BBQueue::size() {
	switch (strategy_) {
	case BranchingStrategy::BFS:
		return bfsQueue_.size();
	case BranchingStrategy::BstFS:
		return bfsPriorityQueue_.size();
	case BranchingStrategy::DFS:
		return dfsStack_.size();
	default:
		return bfsPriorityQueue_.size();
	}
}

bool BBQueue::empty() {
	switch (strategy_) {
	case BranchingStrategy::BFS:
		return bfsQueue_.empty();
	case BranchingStrategy::BstFS:
		return bfsPriorityQueue_.empty();
	case BranchingStrategy::DFS:
		return dfsStack_.empty();
	default:
		return bfsPriorityQueue_.empty();
	}
}

long long int BBQueue::processed() {
	return processed_;
}

void BBQueue::incrementProcessed() {
	++processed_;
}

