/*
 * vector_deque.h
 *
 *  Created on: Jul 20, 2015
 *      Author: breubreubreu
 */

#ifndef VECTOR_DEQUE_H_
#define VECTOR_DEQUE_H_

#include <vector>
#include <fstream>
#include <iostream>
#include <unistd.h>

template <typename T>
class	vector_deque
{
private:
	// Members
	int m_length;

public:
	// Members
	std::vector<T>		m_elements;
	T					dummyLastElement;

	// Constructors
	vector_deque()
	{
		m_length = 0;
	}

	vector_deque(int totalSize)
	{
		m_length = 0;
		m_elements.resize(totalSize);
	}

	// Methods
	void resize(int newLength)
	{
		m_elements.resize(newLength);
	}

	int size()
	{
		return m_length;
	}

	void add(const T& toAdd)
	{
		m_elements[m_length] = toAdd;
		++m_length;
	}

	void remove(int pos)
	{
		if(pos < m_length && m_length > 0)
		{
			dummyLastElement = m_elements[m_length - 1];
			m_elements[pos] = dummyLastElement;
			--m_length;
		}
	}

	void clearAll()
	{
		m_length = 0;
	}

	bool isEmpty()
	{
		return m_length == 0;
	}

    T& operator[](int idx)       { return m_elements[idx]; };
    const T& operator[](int idx) const { return m_elements[idx]; };

    void DataStatus()
    {
    	std::cout << m_length << " : ";
    	for(int iii = 0; iii < m_length; ++iii)
    	{
    		std::cout <<  m_elements[iii]->info().ExtIndex;
    	}
    	std::cout << std::endl;
    }
};

#endif /* VECTOR_DEQUE_H_ */
