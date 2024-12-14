#include <iostream>
#include "zlib.h"
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <queue>
#include <map>
#include <bitset>
#include <cmath>
#include <stdint.h>

using namespace std;

bool compare(const pair<string, float> left, const pair<string, float> right)
{
	return left.second < right.second;
}//функция сортировки

void writeFile(const std::string& filename, const std::string& content) {
	std::ofstream outFile(filename);
	outFile << content; 
	outFile.close();    
}

vector<vector<pair<string, float>>> divid2(vector<pair<string, float>> z) {
	std::sort(begin(z), end(z), compare);
	float halfsum = 0, sum = 0, del2, stop;
	vector<pair<string, float>> u1, u2;
	for (int i = 0; i < z.size(); i++) {
		halfsum += z[i].second;
	}
	halfsum /= 2;
	float del1 = halfsum * 2;
	for (int i = 0; i < z.size(); i++) {
		sum += z[i].second;
		del2 = abs(sum - halfsum);
		if (del2 > del1) {
			stop = i;
			break;
		}
		else {
			u1.push_back(z[i]);
			del1 = del2;
		}
	}
	for (int i = stop; i < z.size(); i++) {
		u2.push_back(z[i]);
	}
	vector<vector<pair<string, float>>> uf = { u1, u2 };
	return uf;

}//функция деления массива на два версия Кузнецова

class block {
public:
	vector<pair<string, float>> u;
	string n;
	block* b1, * b2;
	block(vector<pair<string, float>>& z, string n1) {
		u = z;
		n = n1;
		if (u.size() > 1) {
			vector<vector<pair<string, float>>> uf = divid2(z);
			b1 = new block(uf[0], n + "1");
			b2 = new block(uf[1], n + "0");
		}
		/*/
		else if (u.size() == 1) {
			cout << u[0] << " = " << n << endl;
		}
		//*/
	}
	~block() {
		//delete b1;
		//delete b2;
	}
	map<string, string> getcodes() {
		map<string, string> output;
		vector<string> codes;
		vector<string> z1;
		getcode(codes, z1);
		for (int i = 0; i < z1.size(); i++) {
			output[z1[i]] = codes[i];
		}
		return output;
	}
	void getcode(vector<string>& codes, vector<string>& z1) {
		if (u.size() == 1) {
			codes.push_back(n);
			z1.push_back(u[0].first);
			return;
		}
		b1->getcode(codes, z1);
		b2->getcode(codes, z1);
	}


};//класс дерева для закодирования символов 

struct Node {
	string character;
	float probability;
	Node* left, * right;

	Node(string ch, float prob) : character(ch), probability(prob), left(nullptr), right(nullptr) {}
};


struct Compare1 {
	bool operator()(Node* l, Node* r) {
		return l->probability > r->probability;
	}
};


void generateCodes(Node* root, string code, map<string, string>& codes) {
	if (root == nullptr) return;


	if (!root->left && !root->right) {
		codes[root->character] = code;
	}

	generateCodes(root->left, code + "0", codes);
	generateCodes(root->right, code + "1", codes);
}

string readFile(string name) {
	string line, text = "";

	ifstream in(name);
	in.is_open();
	while (getline(in, line)) {
		text += line;
	}
	in.close();
    std::cout << text << std::endl;
	return text;
}

bool find(vector<pair<string, float>> ensemble, string subtext) {
	for (int i = 0; i < ensemble.size(); i++) {
		if (ensemble[i].first == subtext) {
			return 1;
		}
	}
	return 0;
}

float count(string text, string subtext) {
	float amount = 0;
	for (int i = 0; i < text.size(); i += subtext.size()) {
		if (subtext == text.substr(i, subtext.size())) {
			amount++;
		}
	}
	float fear = amount / (text.size() / subtext.size());
	return fear;
}

vector<pair<string, float>> getEnsemble(string text, int l) {
	vector<pair<string, float>> ensemble;
	string s(l - (text.size() % l), ' '), subtext;
	if (s.size() != l) {
		text += s;
	}

	for (int i = 0; i < text.size(); i += l) {
		subtext = text.substr(i, l);
		if (!find(ensemble, subtext)) {
			ensemble.push_back({ subtext, count(text, subtext) });
		}
	}
	return ensemble;
}

void print(vector<pair<string, float>> ensemble) {
	float fear = 0;
	for (int i = 0; i < ensemble.size(); i++) {
		cout << ensemble[i].first << "   " << ensemble[i].second << endl;
		fear += ensemble[i].second;
	}
	cout << fear << endl;
	float entropy = 0;
	for (int i = 0; i < ensemble.size(); i++) {
		entropy += ensemble[i].second * log(ensemble[i].second) / log(2);
	}
	cout << "Энтропия = " << -entropy << endl;
	cout << "Избыточность = " << 1 + entropy / (log(ensemble.size()) / log(2)) << endl;
	cout << endl;
}

int get_optimal_long_block(string text) {
	float entropy = 0, max_entropi = 0;
	vector<pair<string, float>> optimal_ensemble, ensemble;
	for (int i = 1; i < 11; i++) {
		ensemble = getEnsemble(text, i);
		entropy = 0;
		for (int ii = 0; ii < ensemble.size(); ii++) {
			entropy += ensemble[ii].second * log(ensemble[ii].second) / log(2);
		}
		entropy *= (-1);
		if (entropy > max_entropi) {
			max_entropi = entropy;
		}
		else {
			return (i - 1);
		}
	}
	return 10;
}

void compression_file_1(string name) {
	std::string line = readFile(name);
	int optimal_long_block = get_optimal_long_block(line);
	vector<pair<string, float>> ensemble = getEnsemble(line, optimal_long_block);
	priority_queue<Node*, vector<Node*>, Compare1> minHeap;
	for (int i = 0; i < ensemble.size(); i++) {
		minHeap.push(new Node(ensemble[i].first, ensemble[i].second));
	}
	while (minHeap.size() > 1) {
		Node* left = minHeap.top();
		minHeap.pop();
		Node* right = minHeap.top();
		minHeap.pop();
		Node* newNode = new Node("0", left->probability + right->probability);
		newNode->left = left;
		newNode->right = right;
		minHeap.push(newNode);
	}
	map<string, string> codes;
	generateCodes(minHeap.top(), "", codes);
	int max_cod_long = 0;
	for (const auto& pair : codes) {
		if (pair.second.size() > max_cod_long) {
			max_cod_long = pair.second.size();
		}
	}
	int byts = (max_cod_long / 8) + (bool(max_cod_long % 8));
	for (auto& pair : codes) {
		string cod = "";
		for (int i = 0; i < pair.second.size(); i += 8) {
			int subcod = stoi(pair.second.substr(i, 8), nullptr, 2);
			char ch = char(subcod);
			cod += ch;
		}
		pair.second = cod;
	}
	ofstream out_file{"comp_" + name};
	string s(optimal_long_block - (line.size() % optimal_long_block), ' ');
	if (s.size() != optimal_long_block) {
		line += s;
	}
	for (int i = 0; i < line.size(); i += optimal_long_block) {
		string sublin = line.substr(i, optimal_long_block);
		out_file << codes[sublin];
	}
	out_file << endl << "codes" << byts << " " << optimal_long_block << endl;
	for (const auto& pair : codes) {
		out_file << pair.first << pair.second;
	}
	out_file.close();
}

void compression_file_2(string name) {
	std::string line = readFile(name);
	int optimal_long_block = get_optimal_long_block(line);
	vector<pair<string, float>> ensemble = getEnsemble(line, optimal_long_block);
	priority_queue<Node*, vector<Node*>, Compare1> minHeap;
	for (int i = 0; i < ensemble.size(); i++) {
		minHeap.push(new Node(ensemble[i].first, ensemble[i].second));
	}
	while (minHeap.size() > 1) {
		Node* left = minHeap.top();
		minHeap.pop();
		Node* right = minHeap.top();
		minHeap.pop();
		Node* newNode = new Node("0", left->probability + right->probability);
		newNode->left = left;
		newNode->right = right;
		minHeap.push(newNode);
	}
	map<string, string> codes;
	generateCodes(minHeap.top(), "", codes);
	int min_bit = 2147483647, max_cod_long = 0;
	for (const auto& pair : codes) {
		cout << pair.first << ": " << pair.second << endl;
		if (pair.second.size() < min_bit) {
			min_bit = pair.second.size();
		}
		if (pair.second.size() > max_cod_long) {
			max_cod_long = pair.second.size();
		}
	}
	cout << min_bit << endl;
	string cod, block, bit_cod;
	char ch;
	int j = 7;
	for (int i = 0; i < line.size(); i += optimal_long_block) {
		block = line.substr(i, optimal_long_block);
		cod = codes[block];
		for (int ii = 0; ii < cod.size(); ii++) {
			ch |= ((cod[ii] - 48) << j);
			j--;
			if (j < 0) {
				bit_cod.push_back(ch);
				j = 7;
				ch = 0;
			}
		}
	}
	bit_cod.push_back(ch);
	ofstream out_file{ "comp_" + name };
	out_file << bit_cod << endl;
	out_file << "cod " << min_bit << optimal_long_block << endl;
	int byts = (max_cod_long / 8) + (bool(max_cod_long % 8));
	for (auto& pair : codes) {
		string cod = "";
		for (int i = 0; i < pair.second.size(); i += 8) {
			int subcod = stoi(pair.second.substr(i, 8), nullptr, 2);
			char ch = char(subcod);
			cod += ch;
		}
		pair.second = cod;
	}
	for (const auto& pair : codes) {
		out_file << pair.first << pair.second << " ";
	}
	out_file.close();
}

void compression_file_3(string name) {
	std::string line = readFile(name);
	int optimal_long_block = get_optimal_long_block(line);
	vector<pair<string, float>> ensemble = getEnsemble(line, optimal_long_block);
	block b(ensemble, "");
	map<string, string> codes = b.getcodes();
	cout << codes.size() << endl;
	for (const auto& pair : codes) {
		cout << pair.first << ": " << pair.second << endl;
	}
	int min_bit = 2147483647, max_cod_long = 0;
	for (const auto& pair : codes) {
		cout << pair.first << ": " << pair.second << endl;
		if (pair.second.size() < min_bit) {
			min_bit = pair.second.size();
		}
		if (pair.second.size() > max_cod_long) {
			max_cod_long = pair.second.size();
		}
	}
	cout << min_bit << endl;
	string cod, block, bit_cod;
	char ch;
	int j = 7;
	for (int i = 0; i < line.size(); i += optimal_long_block) {
		block = line.substr(i, optimal_long_block);
		cod = codes[block];
		for (int ii = 0; ii < cod.size(); ii++) {
			ch |= ((cod[ii] - 48) << j);
			j--;
			if (j < 0) {
				bit_cod.push_back(ch);
				j = 7;
				ch = 0;
			}
		}
	}
	bit_cod.push_back(ch);
	ofstream out_file{ "comp_" + name };
	out_file << bit_cod << endl;
	out_file << "cod " << min_bit << optimal_long_block << endl;
	int byts = (max_cod_long / 8) + (bool(max_cod_long % 8));
	for (auto& pair : codes) {
		string cod = "";
		for (int i = 0; i < pair.second.size(); i += 8) {
			int subcod = stoi(pair.second.substr(i, 8), nullptr, 2);
			char ch = char(subcod);
			cod += ch;
		}
		pair.second = cod;
	}
	for (const auto& pair : codes) {
		out_file << pair.first << pair.second << " ";
	}
	out_file.close();
}

string replacement(map<string, string> &codes, string &text, int &optimal_long_block, int max_cod_long) {
	string s(optimal_long_block - (text.size() % optimal_long_block), ' '), block, output = "", cod;
	if (s.size() != optimal_long_block) {
		text += s;
	}
	for (int i = 0; i < text.size(); i += optimal_long_block) {
		block = text.substr(i, optimal_long_block);
		output += codes[block];
	}
	if (output.length() % 8 != 0) {
		string padding(8 - (output.length() % 8), '0');
		output += padding;
	}
	bitset<32> b1(output.size());
	string sizetext = b1.to_string();
	bitset<32> b2(optimal_long_block);
	sizetext += b2.to_string();
	bitset<32> b3(max_cod_long);
	sizetext += b3.to_string();
	output = sizetext + output;
	
	for (const auto& pair : codes) {
		block = pair.first;
		cod = pair.second;
		int long_cod = cod.size();
		bitset<8> b3(long_cod);
		for (char c : block) {
			output += (bitset<8>(c).to_string());
		}
		output += b3.to_string();
		output += cod;
	}

	return output;
}

void writeBinaryToFile(std::string binaryString, const std::string& fileName) {
	std::ofstream outFile(fileName, std::ios::binary);
	std::string::size_type len = binaryString.length();
	for (std::string::size_type i = 0; i < len; i += 8) {
		std::bitset<8> bits(binaryString.substr(i, 8));
		unsigned char byte = static_cast<unsigned char>(bits.to_ulong());
		outFile.write(reinterpret_cast<const char*>(&byte), sizeof(byte));
	}
	outFile.close();
}

int long_block_Shannon(string text) {
	int long_block = get_optimal_long_block(text), min_size_text = INT_MAX, best_length;
	for (int i = 1; i <= long_block; i++) {
		int optimal_long_block = i;
		vector<pair<string, float>> ensemble = getEnsemble(text, optimal_long_block);
		block b(ensemble, "");
		map<string, string> codes = b.getcodes();
		string cod_text = replacement(codes, text, optimal_long_block, 1);
		if (cod_text.size() < min_size_text) {
			min_size_text = cod_text.size();
			best_length = i;

		}
		cout << i << " " << cod_text.size() << endl;
	}
	return best_length;
}

void compression_file_Shannon(string name) {
	string line = readFile(name);
	int optimal_long_block = long_block_Shannon(line);
	//int optimal_long_block = 1;
	vector<pair<string, float>> ensemble = getEnsemble(line, optimal_long_block);
	block b(ensemble, "");
	map<string, string> codes = b.getcodes();
	cout << codes.size() << endl;
	int min_bit = 2147483647, max_cod_long = 0;
	for (const auto& pair : codes) {
		cout << pair.first << ": " << pair.second << endl;
		if (pair.second.size() < min_bit) {
			min_bit = pair.second.size();
		}
		if (pair.second.size() > max_cod_long) {
			max_cod_long = pair.second.size();
		}
	}
	cout << max_cod_long << endl;
	if (max_cod_long % 8 != 0) {
		max_cod_long += 8 - (max_cod_long % 8);
	}
	max_cod_long /= 8;
	cout << max_cod_long << endl;
	string cod_text = replacement(codes, line, optimal_long_block, max_cod_long);
	string new_name = name.substr(0, name.size() - 4) + ".dat";
	cout << new_name << endl;
	writeBinaryToFile(cod_text, new_name);
}

string readBinaryFromFile(const std::string& fileName) {
	std::ifstream inFile(fileName, std::ios::binary);
	std::string binaryString;
	unsigned char byte;
	while (inFile.read(reinterpret_cast<char*>(&byte), sizeof(byte))) {
		std::bitset<8> bits(byte);
		binaryString += bits.to_string();
	}
	inFile.close();
	return binaryString;
}

void decompression(string name) {
	string input = readBinaryFromFile(name);
	int size_text = std::stoi(input.substr(0, 32), nullptr, 2);
	int long_block = std::stoi(input.substr(32, 32), nullptr, 2);
	int byt_cod = std::stoi(input.substr(64, 32), nullptr, 2);
	cout << size_text << " " << long_block << " " << byt_cod << endl;
	string s_codes = input.substr(size_text + 96, input.size() - size_text - 96);
	map<string, string> codes;
	for (int i = 0; i < s_codes.size();) {
		string bit_block = s_codes.substr(i, long_block * 8);
		int long_cod = stoi(s_codes.substr(i + long_block * 8, 8), nullptr, 2);
		std::string block;
		for (size_t ii = 0; ii < bit_block.length(); ii += 8) {
			string byte = bit_block.substr(ii, 8);
			char character = static_cast<char>(std::stoi(byte, nullptr, 2));
			block += character;
		}
		string cod = s_codes.substr(i + long_block * 8 + 8, long_cod);
		codes[cod] = block;
		i += long_block * 8 + 8 + long_cod;
	}
	for (const auto& pair : codes) {
		cout << pair.first << ": " << pair.second << endl;
	}
	string cod_text = input.substr(96, size_text);
	string text = "", podtext = "";
	for (int i = 0; i < cod_text.size(); i++) {
		podtext.push_back(cod_text[i]);
		if (codes.count(podtext)) {
			text += codes[podtext];
			podtext = "";
		}
	}
	cout << text << endl;
	string newname = name.substr(0, name.size() - 4) + "1.txt";
	writeFile(newname, text);
}

void compression_file_Huffman(string name) {
	std::string line = readFile(name);
	int optimal_long_block = long_block_Shannon(line);

	vector<pair<string, float>> ensemble = getEnsemble(line, optimal_long_block);
	priority_queue<Node*, vector<Node*>, Compare1> minHeap;
	for (int i = 0; i < ensemble.size(); i++) {
		minHeap.push(new Node(ensemble[i].first, ensemble[i].second));
	}
	while (minHeap.size() > 1) {
		Node* left = minHeap.top();
		minHeap.pop();
		Node* right = minHeap.top();
		minHeap.pop();
		Node* newNode = new Node("0", left->probability + right->probability);
		newNode->left = left;
		newNode->right = right;
		minHeap.push(newNode);
	}
	map<string, string> codes;
	generateCodes(minHeap.top(), "", codes);
	int min_bit = 2147483647, max_cod_long = 0;
	for (const auto& pair : codes) {
		cout << pair.first << ": " << pair.second << endl;
		if (pair.second.size() < min_bit) {
			min_bit = pair.second.size();
		}
		if (pair.second.size() > max_cod_long) {
			max_cod_long = pair.second.size();
		}
	}
	if (max_cod_long % 8 != 0) {
		max_cod_long += 8 - (max_cod_long % 8);
	}
	max_cod_long /= 8;
	cout << max_cod_long << endl;
	string cod_text = replacement(codes, line, optimal_long_block, max_cod_long);
	string new_name = name.substr(0, name.size() - 4) + ".dat";
	cout << new_name << endl;
	writeBinaryToFile(cod_text, new_name);
}

int main() {
	setlocale(LC_ALL, "RU");
	/*/
	std::string line = readFile();
	std::cout << line << std::endl;
	vector<pair<string, float>> ensemble = getEnsemble(line, 1);
	print(ensemble);
	vector<pair<string, float>> ensemble1 = getEnsemble(line, 2);
	print(ensemble1);
	vector<pair<string, float>> ensemble2 = getEnsemble(line, 3);
	print(ensemble2);
	vector<pair<string, float>> ensemble3 = getEnsemble(line, 4);
	print(ensemble3);
	vector<pair<string, float>> ensemble4 = getEnsemble(line, 5);
	print(ensemble4);

	//*/
	/*/
	std::string line = readFile("text.txt");
	std::cout << line << std::endl;
	//int optimal_long_block = get_optimal_long_block(line);
	int optimal_long_block = 2;
	cout << optimal_long_block << endl;
	vector<pair<string, float>> ensemble = getEnsemble(line, optimal_long_block);
	print(ensemble);
	priority_queue<Node*, vector<Node*>, Compare1> minHeap;
	for (int i = 0; i < ensemble.size(); i++) {
		minHeap.push(new Node(ensemble[i].first, ensemble[i].second));
	}
	while (minHeap.size() > 1) {
		Node* left = minHeap.top();
		minHeap.pop();
		Node* right = minHeap.top();
		minHeap.pop();

		Node* newNode = new Node("0", left->probability + right->probability);
		newNode->left = left;
		newNode->right = right;
		minHeap.push(newNode);
	}
	map<string, string> codes;
	generateCodes(minHeap.top(), "", codes);
	cout << codes.size() << endl;
	float length = 0;
	int j = 0, max_cod_long = 0;
	for (const auto& pair : codes) {
		cout << pair.first << ": " << pair.second << endl;
		length += ensemble[j].second * pair.second.size();
		if (pair.second.size() > max_cod_long) {
			max_cod_long = pair.second.size();
		}
		j++;
	}
	cout << length << endl;
	cout << max_cod_long << endl;
	int byts = (max_cod_long / 8) + (bool(max_cod_long % 8));
	cout << byts << endl;
	for (auto& pair : codes) {
		string cod = "";
		for (int i = 0; i < pair.second.size(); i += 8) {
			int subcod = stoi(pair.second.substr(i, 8), nullptr, 2);
			char ch = char(subcod);
			cod += ch;
		}
		pair.second = cod;
	}
	for (const auto& pair : codes) {
		cout << pair.first << ": " << pair.second << endl;
	}
	ofstream out_file{ "copy.txt" };
	string s(optimal_long_block - (line.size() % optimal_long_block), ' ');
	if (s.size() != optimal_long_block) {
		line += s;
	}
	for (int i = 0; i < line.size(); i += optimal_long_block) {
		string sublin = line.substr(i, optimal_long_block);
		out_file << codes[sublin];
	}
	out_file << endl << "codes" << byts << " " << optimal_long_block << endl;
	for (const auto& pair : codes) {
		out_file << pair.first << pair.second;
	}
	out_file.close();
	//*/

	cout << "Введите имя файла" << endl;
	string name;
	cin >> name;
	cout << "Вы хотите заархивировать(1) файл, или разорхевировать(2)?" << endl;
	int answer;
	cin >> answer;
	if (answer == 1) {
		cout << "Вы хотите заархивировать с помощью алгоритмом Шеннона - Фано(1) или Хаффмана(2)" << endl;
		cin >> answer;
		if (answer == 1) {
			compression_file_Shannon(name);
		}
		else {
			compression_file_Huffman(name);
		}
	}
	else {
		decompression(name);
	}
	return 0;
}
