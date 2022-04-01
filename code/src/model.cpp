#include <boost/functional/hash.hpp>
#include <boost/bind/bind.hpp>
#include <iostream>
#include <fstream>
#include <list>
#include <stdexcept>
#include <random>
#include <algorithm>
#include <chrono>

#include "model.hpp"

namespace dotl {

  bool in(const Proposition& p, const Propositions& P) {
    return P.contains(p);
  };

  bool Valuate::operator() (const Proposition& p) const {
    return in(p, (this->D).get_true_propositions(this->h));
  };
  bool Valuate::operator() (const Top&) const { return true; };
  bool Valuate::operator() (const Bot&) const { return false; };
  bool Valuate::operator() (const Neg& phi) const {
    return !boost::apply_visitor(*this, *phi.phi);
  };
  bool Valuate::operator() (const Con& phi) const {
    return boost::apply_visitor(*this, *phi.phi) &&
           boost::apply_visitor(*this, *phi.psi);
  };
  bool Valuate::operator() (const Dis& phi) const {
    return boost::apply_visitor(*this, *phi.phi) ||
           boost::apply_visitor(*this, *phi.psi);
  };
  bool Valuate::operator() (const Believes& phi) const {
    auto& D = this->D;
    auto& h = this->h;
    auto visitor = Valuate(D, h);
    for (auto& h2 : D.get_minimal_histories(phi.i, h)) {
      visitor.h = h2;
      if (D.has_edge(phi.i, h, h2) && !boost::apply_visitor(visitor, *phi.phi))
      {
          return false;
      };
    };
    return true;
  };
  bool Valuate::operator() (const Knows& phi) const {
    auto& D = this->D;
    auto& h = this->h;
    auto visitor = Valuate(D, h);
    for (auto& h2 : D.edges_from(phi.i, h)) {
      visitor.h = h2;
      if (!boost::apply_visitor(visitor, *phi.phi))
      {
          return false;
      };
    };
    for (auto& h2 : D.edges_to(phi.i, h)) {
      visitor.h = h2;
      if (!boost::apply_visitor(visitor, *phi.phi))
      {
          return false;
      };
    };
    return true;
  };
  bool Valuate::operator() (const Common_knowledge& phi) const {
    auto& D = this->D;
    auto& h = this->h;
    auto visitor = Valuate(D, h);
    for (auto& h2 : D.get_group_reachable(h)) {
      visitor.h = h2;
      if (!boost::apply_visitor(visitor, *phi.phi))
      {
          return false;
      };
    };
    return true;
  };
  bool Valuate::operator() (const After_event& phi) const {
    return boost::apply_visitor(With_path(*this, *phi.phi), *phi.pi);
  };

  bool With_path::operator() (const Event_id& e) const {
    auto& D = this->D;
    auto& child = D.get_child(this->h, e);
    if (child.has_value())
      return D.valuate(this->phi, child.value());
    return true;
  };
  bool With_path::operator() (const Seq& pi) const {
    return this->D.valuate(
        After_event(After_event(this->phi, *pi.pi2), *pi.pi1), 
        this->h);
  };
  bool With_path::operator() (const Or& pi) const {
    auto& D = this->D;
    auto& h = this->h;
    return D.valuate(After_event(this->phi, *pi.pi1), h) ||
      D.valuate(After_event(this->phi, *pi.pi2), h);
  };
  bool With_path::operator() (const Inv& pi) const {
    return boost::apply_visitor(
        With_inv_path(this->D, this->h, this->phi), *pi.pi);
  };

  bool With_inv_path::operator() (const Event_id& e) const {
    auto& D = this->D;
    auto& parent = D.get_parent(this->h, e);
    if (parent.has_value())
      return D.valuate(this->phi, parent.value());
    return false;
  };
  bool With_inv_path::operator() (const Seq& pi) const {
    return (this->D).valuate(
        After_event(After_event(this->phi, Inv(*pi.pi1)), Inv(*pi.pi2)),
        this->h
        );
  };
  bool With_inv_path::operator() (const Or& pi) const {
    auto& D = this->D;
    auto& h = this->h;
    return D.valuate(After_event(this->phi, Inv(*pi.pi2)), h) ||
      D.valuate(After_event(this->phi, Inv(*pi.pi1)), h);
  };
  bool With_inv_path::operator() (const Inv& pi) const {
    return boost::apply_visitor(
        With_path(this->D, this->h, this->phi), *pi.pi);
  };

  std::string Formula_to_string::operator() (const Proposition& p) const {
    std::string s;
    s = p;
    return s;
  };
  std::string Formula_to_string::operator() (const Top&) const {
    return "T";
  };
  std::string Formula_to_string::operator() (const Bot&) const {
    return "F";
  };
  std::string Formula_to_string::operator() (const Neg& phi) const {
    auto s = boost::apply_visitor(Formula_to_string(), *phi.phi);
    return "~" + s;
  };
  std::string Formula_to_string::operator() (const Con& phi) const {
    auto s1 = boost::apply_visitor(Formula_to_string(), *phi.phi);
    auto s2 = boost::apply_visitor(Formula_to_string(), *phi.psi);
    return s1 + " /\\ " + s2;
  };
  std::string Formula_to_string::operator() (const Dis& phi) const {
    auto s1 = boost::apply_visitor(Formula_to_string(), *phi.phi);
    auto s2 = boost::apply_visitor(Formula_to_string(), *phi.psi);
    return s1 + " \\/ " + s2;
  };
  std::string Formula_to_string::operator() (const Believes& phi) const {
    auto s = boost::apply_visitor(Formula_to_string(), *phi.phi);
    return "B" + std::to_string(phi.i) + " " + s;
  };
  std::string Formula_to_string::operator() (const Knows& phi) const {
    auto s = boost::apply_visitor(Formula_to_string(), *phi.phi);
    return "K" + std::to_string(phi.i) + " " + s;
  };
  std::string Formula_to_string::operator() (const Common_knowledge& phi) const {
    auto s = boost::apply_visitor(Formula_to_string(), *phi.phi);
    return "C " + s;
  };
  std::string Formula_to_string::operator() (const After_event& phi) const {
    auto s1 = boost::apply_visitor(Path_to_string(), *phi.pi);
    auto s2 = boost::apply_visitor(Formula_to_string(), *phi.phi);
    return "<" + s1 + ">" + s2;
  };
  std::string Path_to_string::operator() (const Event_id& e) const {
    return "e" + std::to_string(e);
  };
  std::string Path_to_string::operator() (const Seq& pi) const {
    auto s1 = boost::apply_visitor(Path_to_string(), *pi.pi1);
    auto s2 = boost::apply_visitor(Path_to_string(), *pi.pi2);
    return "(" + s1 + ";" + s2 + ")";
  };
  std::string Path_to_string::operator() (const Or& pi) const {
    auto s1 = boost::apply_visitor(Path_to_string(), *pi.pi1);
    auto s2 = boost::apply_visitor(Path_to_string(), *pi.pi2);
    return "(" + s1 + "/\\" + s2 + ")";
  };
  std::string Path_to_string::operator() (const Inv& pi) const {
    auto s = boost::apply_visitor(Path_to_string(), *pi.pi);
    return "(" + s + ")" + "^-1";
  };
  std::ostream& operator<<(std::ostream& os, const Formula& phi) {
    return os << boost::apply_visitor(Formula_to_string(), phi);
  };
  std::ostream& operator<<(std::ostream& os, const Path& pi) {
    return os << boost::apply_visitor(Path_to_string(), pi);
  };

  std::vector<Id>& edges(Relation& R, Agent_id i, Id from) {
    return R[ {i, from} ].second;
  };

  const bool has_edge(Relation& R, Agent_id i, Id from, Id to) {
    auto& [sorted, one_reachable] = R.at( {i, from} );
    if (!sorted) {
      std::sort(one_reachable.begin(), one_reachable.end());
      sorted = true;
    };
    return std::binary_search(one_reachable.begin(), one_reachable.end(), to);
  };

  Event_id Event_model::total_events = 0;

  const Propositions& Event_model::get_positive_effects(Event_id event_id) const {
    return this->events[event_id-this->offset].positive_effects;
  };

  const Propositions& Event_model::get_negative_effects(Event_id event_id) const {
    return this->events[event_id-this->offset].negative_effects;
  };

  const Formula& Event_model::get_precondition(Event_id event_id) const {
    return this->events[event_id-this->offset].precondition;
  };

  bool Event_model::is_designated(Event_id event_id) const {
    return this->events[event_id-this->offset].designated;
  }

  const bool Event_model::has_edge(Agent_id i, Event_id from, Event_id to) {
    return dotl::has_edge(this->plausibility_relation, i, from, to);
  };
  const bool Event_model::indistinguishable_from(Agent_id i, Event_id e1,
      Event_id e2) {
    return this->has_edge(i, e1, e2) || this->has_edge(i, e2, e1);
  };

  void Event_model::perspective_shift(Agent_id agent_id) {
    std::list<Event_id> designated;
    for (auto e : this->events) {
      if (e.designated) {
        designated.push_back(e.id);
      };
    };
    Event_id from;
    while (!designated.empty()) {
      from = designated.front();
      designated.pop_front();
      for (auto& to : this->events) {
        if (!to.designated && this->indistinguishable_from(agent_id, from, to.id)) {
          to.designated = true;
          designated.push_back(to.id);
        };
      };
    };
  };

  void add_edge_to_relation(Relation& relation, Agent_id agent_id, Edge edge) {
    relation[ {agent_id, edge.first} ].second.emplace_back(edge.second);
  };

  Event_id Event_model::new_event(bool designated, 
      Formula& precondition,
      std::initializer_list<Proposition> positive_effects,
      std::initializer_list<Proposition> negative_effects) {
    auto id = total_events++;
    events.emplace_back(Event(id, designated, precondition, positive_effects,
          negative_effects));
    for (int i = 0; i < this->agents; i++) {
      add_edge_to_relation(this->plausibility_relation, i, {id, id});
      this->union_relation.insert({id,id});
    };
    return id;
  };

  void add_edges_to_relation(Relation& relation, Agent_id agent_id,
      std::vector<Edge> edges) {
    for (Edge edge : edges) add_edge_to_relation(relation, agent_id, edge);
  };
  void add_inverse_edges_to_relation(Relation& relation, Agent_id agent_id,
      std::vector<Edge> edges) {
    for (Edge edge : edges) 
      add_edge_to_relation(relation, agent_id, {edge.second, edge.first});
  };

  void Event_model::add_plausibility_relation(Agent_id agent_id,
      std::vector<Edge> edges) {
    add_edges_to_relation(this->plausibility_relation, agent_id, edges);
    // Maintain transitive closure.
    auto& tc = this->union_relation;
    for (auto& [a,b] : edges) {
      if (!tc.contains({a,b})) {
        tc.insert({a,b});
        for (auto& e1 : this->events) {
          for (auto& e2 : this->events) {
            if (tc.contains({e1.id,a}) && tc.contains({b,e2.id})) {
              tc.insert({e1.id,e2.id});
            };
          };
        };
      };
    };
  };

  void Event_model::print() {
    for (const Event& event : events) {
      std::cout << "e" << event.id << " : <";
      std::cout << event.precondition;
      std::cout << ", ";
      auto first = true;
      for (const Proposition& p : event.positive_effects) {
        if (first) {
          std::cout << p;
          first = false;
        }
        else std::cout << " /\\" << p;
      }
      first = true;
      for (const Proposition& p : event.negative_effects) {
        if (first) {
          std::cout << "~" << p;
          first = false;
        }
        else std::cout << " /\\ ~" << p;
      }
      if (event.negative_effects.empty() && event.positive_effects.empty()) {
        std::cout << "T";
      };
      std::cout << ">" << std::endl;
      for (auto agent = 0; agent < this->agents; agent++) {
        std::cout << "<=_" << agent << ": {";
        bool separator = false;
        for (const Event_id most_plausible : 
            this->plausibility_relation[ {agent, event.id} ].second) {
          if (separator) { 
            std::cout << ",";
          }
          else {
            separator = true;
          };
          std::cout << "(e" << event.id << "," << "e" << most_plausible << ")";
        }
        std::cout << "}" << std::endl;
      }
    }
  };

  bool History::operator<(const History& other) const {
    return this->id < other.id;
  };

  const Propositions& DoTL_model::get_true_propositions(History_id history_id) const {
    return this->histories[history_id].true_propositions;
  };

  std::vector<History_id> DoTL_model::get_group_reachable(
      History_id history_id) {
    std::vector<History_id> reachable;
    auto length = this->histories[history_id].length;
    std::unordered_set<History_id> visited;
    auto offset = histories_by_length[length][0];
    visited.emplace(history_id);
    std::list<size_t> queue;
    queue.push_back(history_id);

    while (!queue.empty()) {
      history_id = queue.front();
      queue.pop_front();
      reachable.emplace_back(history_id);

      for (auto i = 0; i < this->agents; i++) {
        std::vector<History_id>& relation = this->edges_from(i, history_id);
        for (auto it = relation.begin(); it != relation.end(); ++it) {
          if (!visited.contains(*it)) {
            visited.emplace(*it);
            queue.push_back(*it);
          };
        };
        std::vector<History_id>& inv_relation = this->edges_to(i, history_id);
        for (auto it = inv_relation.begin(); it != inv_relation.end(); ++it) {
          if (!visited.contains(*it)) {
            visited.emplace(*it);
            queue.push_back(*it);
          };
        };
      };
    };

    return reachable;
  };

 std::vector<History_id> DoTL_model::get_information_cell(Agent_id agent_id,
     History_id history_id) {
    std::vector<History_id> information_cell;
    auto length = this->histories[history_id].length;
    std::unordered_set<History_id> visited;
    auto offset = histories_by_length[length][0];
    visited.emplace(history_id);
    std::list<size_t> queue;
    queue.push_back(history_id);

    while (!queue.empty()) {
      history_id = queue.front();
      queue.pop_front();
      information_cell.emplace_back(history_id);

      std::vector<History_id>& relation = this->edges_from(agent_id,
          history_id);
      for (auto it = relation.begin(); it != relation.end(); ++it) {
        if (!visited.contains(*it)) {
          visited.emplace(*it);
          queue.push_back(*it);
        };
      };
      std::vector<History_id>& inv_relation = this->edges_to(agent_id,
          history_id);
      for (auto it = inv_relation.begin(); it != inv_relation.end(); ++it) {
        if (!visited.contains(*it)) {
          visited.emplace(*it);
          queue.push_back(*it);
        };
      };
    };

    return information_cell;
  };

  std::vector<History_id> DoTL_model::get_minimal_histories(Agent_id agent_id,
      std::vector<History_id>& information_cell) {
    std::vector<History_id> minimal;
    for (auto& h : information_cell) {
      bool is_minimal = true;
      int i = 0;
      while (i < information_cell.size() && is_minimal) {
        is_minimal =
          is_minimal && this->has_edge(agent_id, information_cell[i++], h);
      };
      if (is_minimal) minimal.emplace_back(h);
    };
    return minimal;
  };

  std::vector<History_id> DoTL_model::get_minimal_histories(Agent_id agent_id,
      History_id history_id) {
    auto information_cell = this->get_information_cell(agent_id, history_id);
    return this->get_minimal_histories(agent_id, information_cell);
  };

  std::vector<History_id>& DoTL_model::edges_from(Agent_id
      agent_id, History_id history_id) {
    return
      dotl::edges(this->plausibility_relation, agent_id, history_id);
  };
  std::vector<History_id>& DoTL_model::edges_to(Agent_id
      agent_id, History_id history_id) {
    return
      dotl::edges(this->inverse_plausibility_relation, agent_id, history_id);
  };

  const std::optional<History_id> DoTL_model::get_child(History_id history_id,
      Event_id event_id) const {
    auto child_length = this->histories[history_id].length+1;
    if (this->histories_by_length.size() > child_length)
      for (auto child_id : histories_by_length[child_length]) {
        auto& child = this->histories[child_id];
        if (child.parent == history_id && child.event == event_id)
          return child_id;
      };
    return {};
  };

  const std::optional<History_id> DoTL_model::get_parent(History_id history_id,
      Event_id event_id) const {
    if (this->histories.size() > history_id && 
        this->histories[history_id].event == event_id)
      return this->histories[history_id].parent;
    return {};
  };

  bool DoTL_model::is_designated(History_id history_id) const {
    return this->histories[history_id].designated;
  }

  bool DoTL_model::valuate(const Formula& phi, History_id history_id) {
    return boost::apply_visitor(Valuate(*this, history_id), phi);
  };

  bool DoTL_model::valuate(const Formula& phi) {
    bool result = true;
    for (const History& h : this->histories_by_length[max_length]) {
      if (h.designated) result = result && this->valuate(phi, h.id);
    }
    return result;
  };

  const bool DoTL_model::has_edge(Agent_id i, History_id from, History_id to) {
      return dotl::has_edge(this->plausibility_relation, i, from, to);
  };
  const bool DoTL_model::indistinguishable_from(Agent_id i, History_id h1,
      History_id h2) {
    return this->has_edge(i, h1, h2) || this->has_edge(i, h2, h1);
  };

  void DoTL_model::add_plausibility_relation(Agent_id agent_id,
      std::vector<Edge> edges) {
    add_edges_to_relation(this->plausibility_relation, agent_id, edges);
    add_inverse_edges_to_relation(
        this->inverse_plausibility_relation, agent_id, edges);
    };

  void DoTL_model::add_history(const History& h) {
    this->histories.emplace_back(h);
    if (this->histories_by_length.size() <= h.length)
      this->histories_by_length.emplace_back(std::vector<History_id>());
    this->histories_by_length[h.length].emplace_back(h.id);
  };

  DoTL_model::DoTL_model(
      std::size_t agents, std::initializer_list<Proposition> true_propositions)
  {
    this->agents = agents;
    this->max_length = 0;
    auto id = this->histories.size();
    this->add_history(History(id, true, true_propositions, id, 0));
    for (int i = 0; i < this->agents; i++) {
      this->add_plausibility_relation(i, { {id, id} });
    };
    this->union_relation.insert( {id, id} );
    this->rng.seed(666);
  };

  bool DoTL_model::find_history(std::vector<History>& hs, History& h) {
    return std::binary_search(hs.begin(), hs.end(), h);
  };
  bool DoTL_model::find_history_id(std::vector<History_id>& is, History_id& i)
  {
    return std::binary_search(is.begin(), is.end(), i);
  };

  std::vector<History_id>& DoTL_model::frontier() {
    return this->histories_by_length[this->max_length];
  };

  // Assumes I is sorted.
  bool DoTL_model::restricted_update(Event_model& event_model,
      std::vector<History_id>* I = NULL) {
    this->sequence_of_events.emplace_back(event_model);
    bool added_history = false;
    bool added_designated_history = false;

    auto first_id = this->histories.size();
    for (History_id hi : this->frontier()) {
      // add_history potentially changes memory location of this->histories.
      History h = this->histories[hi];
      if (I == NULL || this->find_history_id(*I, h.id)) {
        for (auto &e : event_model.events) {
          if (this->valuate(e.precondition, h.id)) {
            // (h,e) and V'(h,e), possible add (h,e) to H_d, and add
            // transition to parent and event that brought on the situation.
            auto id = this->histories.size();
            bool designated = h.designated && e.designated;
            auto he = 
              History(id, designated, h.true_propositions, h.id, h.length+1);
            he.event = e.id;
            for (const Proposition& p : e.positive_effects)
              he.true_propositions.insert(p);
            for (const Proposition& p : e.negative_effects)
              he.true_propositions.erase(p);
            add_history(he);
            added_history = true;
            added_designated_history = 
              added_designated_history || designated;
          };
        };
      }
      else {
        this->marked.emplace_back(h.id);
      };
    };
    // compute plausibility relations.
    int last_id = this->histories.size();
    for (int agent = 0; agent < this->agents; agent++) {
      for (size_t i = first_id; i < last_id; i++) {
        this->add_plausibility_relation(agent, { {i, i} });
        for (size_t j = i+1; j < last_id; j++) {
          History_id hi = this->histories[i].parent;
          History_id hj = this->histories[j].parent;
          Event_id ei = this->histories[i].event;
          Event_id ej = this->histories[j].event;
          bool hihj = 
            this->has_edge(agent, hi, hj);
          bool hjhi =
            this->has_edge(agent, hj, hi);
          bool eiej =
            event_model.has_edge(agent, ei, ej);
          bool ejei =
            event_model.has_edge(agent, ej, ei);
          if (((hihj || hjhi) && eiej && !ejei) ||
              (eiej && ejei && hihj)) {
            this->add_plausibility_relation(agent, { {i, j} });
          };
          if (((hihj || hjhi) && !eiej && ejei) ||
              (eiej && ejei && hjhi)) {
            this->add_plausibility_relation(agent, { {j, i} });
          };
        };
      };
    };
    // compute union of plausibility relations.
    for (size_t i = first_id; i < last_id; i++) {
      this->union_relation.insert( {i,i} );
      for (size_t j = i+1; j < last_id; j++) {
        History_id hi = this->histories[i].parent;
        History_id hj = this->histories[j].parent;
        Event_id ei = this->histories[i].event;
        Event_id ej = this->histories[j].event;
        bool hihj = this->union_relation.contains({hi, hj});
        bool hjhi = this->union_relation.contains({hj, hi});
        bool eiej = event_model.union_relation.contains({ei, ej});
        bool ejei = event_model.union_relation.contains({ej, ei});
        if (((hihj || hjhi) && eiej && !ejei) ||
            (eiej && ejei && hihj)) {
          this->union_relation.insert({i, j});
        };
        if (((hihj || hjhi) && !eiej && ejei) ||
            (eiej && ejei && hjhi)) {
          this->union_relation.insert({j, i});
        };
      };
    };

    if (added_history) this->max_length++;
    return added_designated_history;
  };

  void DoTL_model::action_priority_update(Event_model& event_model) {
    this->restricted_update(event_model);
  };

  // Follows Algorithm 1 from section 3.1.2.
  bool DoTL_model::most_plausible_update(
      Event_model& event_model, Agent_id agent_id) {

    // compute min(agent_id).
    auto& Hmax = this->histories_by_length[this->max_length];
    std::vector<History_id> min_agent_id;
    int num_designated = 0;
    History_id hd;
    for (auto& h : Hmax) {
      if (this->is_designated(h)) {
        num_designated++;
        hd = h;
      };
    };
    if (num_designated) {
      auto hdi = get_information_cell(agent_id, hd);
      if (hdi.size() == num_designated)
        // get_minimal_histories works for general DoTL models. If this proves
        // to be a bottleneck we should optimize get_minimal_histories under
        // the assumption that the model is subjective.
        min_agent_id = get_minimal_histories(agent_id, hdi);
    };

    // compute close(min(agent_id))
    auto& tc = this->union_relation;
    std::vector<History_id> I;
    for (auto h1 : Hmax) {
      for (auto& h2 : min_agent_id) {
        if (tc.contains({ h2, h1 })) {
          I.emplace_back(h1);
          break;
        };
      };
    };
    std::sort(I.begin(), I.end());
    return this->restricted_update(event_model, &I);
  };

  // Assumes marked is sorted by length.
  // Follows Algorithm 4 from section 3.1.2.
  std::list<History_id>::iterator DoTL_model::most_plaus_marked_from_sorted(
      Agent_id agent_id) {
    auto len = [&](auto& h) { return this->histories[h].length; };
    auto designated = [&](auto& h) { return this->is_designated(h); };
    auto begin = std::find_if(this->marked.begin(), this->marked.end(), designated);
    // lambda captures reference to begin, so should get up-to-date version of
    // begin when invoked.
    auto neq_begin_length = [&](auto& h) { return len(h) != len(*begin); };
    while (begin != this->marked.end()) {
      auto end = std::find_if(begin, this->marked.end(), neq_begin_length);
      auto hbegin = this->get_information_cell(agent_id, *begin);
      auto n_designated = std::count_if(
          this->histories_by_length[len(*begin)].begin(),
          this->histories_by_length[len(*begin)].end(),
          designated);
      if (hbegin.size() == n_designated) { // There is a most plausible history.
        auto designated_implies_leq = [&](auto& to, auto& from) {
          return !this->is_designated(from) || 
            this->has_edge(agent_id, from, to);
        };
        auto minimal = [&](auto& h) {
          return std::all_of(begin, end, 
              std::bind_front(designated_implies_leq, h));
        };
        auto h = std::find_if(begin, end, minimal);
        if (h != this->marked.end()) return h;
      };
      begin = std::find_if(end, this->marked.end(), designated);
    };
    throw "Abduction error: no most plausible history amongst marked \
      histories.";
  };

  // The following two abduction heuristics follows Algorithm 3 from section
  // 3.1.2. To increase performance, the implementation of uniform expansion
  // partitions marked histories in place by sorting and searches partitions
  // uniformly at random. All three assume there is at least one marked
  // history.
  std::list<History_id>::iterator DoTL_model::chronological_minimization(
      Agent_id agent_id) {
    auto len = [&](auto& h) { return this->histories[h].length; };
    auto longer = [&](auto& h1, auto& h2) { return len(h1) > len(h2); };
    this->marked.sort(longer);
    return this->most_plaus_marked_from_sorted(agent_id);
  };

  std::list<History_id>::iterator DoTL_model::inv_chronological_minimization(
      Agent_id agent_id) {
    auto len = [&](auto& h) { return this->histories[h].length; };
    auto shorter = [&](auto& h1, auto& h2) { return len(h1) < len(h2); };
    this->marked.sort(shorter);
    return this->most_plaus_marked_from_sorted(agent_id);
  };

  std::list<History_id>::iterator DoTL_model::uniform_expansion(
      Agent_id agent_id) {
    auto len = [&](auto& h) { return this->histories[h].length; };
    auto shorter = [&](auto& h1, auto& h2) { return len(h1) < len(h2); };
    this->marked.sort(shorter);

    std::vector<std::pair<
      std::list<History_id>::iterator, std::list<History_id>::iterator>>
        ranges;
    auto neq_length = [&](auto& h1, auto& h2) { return len(h1) != len(h2); };
    auto begin = this->marked.begin();
    while (begin != this->marked.end()) {
      auto end = std::adjacent_find(begin, this->marked.end(), neq_length);
      if (end != this->marked.end())
        ++end;
      ranges.push_back({begin, end});
      begin = end;
    };
    std::shuffle(ranges.begin(), ranges.end(), rng);
    auto designated = [&](auto& h) { return this->is_designated(h); };
    for (auto& [begin, end] : ranges) {
      begin = std::find_if(begin, end, designated);
      auto hbegin = this->get_information_cell(agent_id, *begin);
      auto n_designated = std::count_if(
          this->histories_by_length[len(*begin)].begin(),
          this->histories_by_length[len(*begin)].end(),
          designated);
      if (hbegin.size() == n_designated) { // There is a most plausible history.
        auto designated_implies_leq = [&](auto& to, auto& from) {
          return !this->is_designated(from) || 
            this->has_edge(agent_id, from, to);
        };
        auto minimal = [&, begin = begin, end = end](auto& h) {
          return std::all_of(begin, end, 
              std::bind_front(designated_implies_leq, h));
        };
        auto h = std::find_if(begin, end, minimal);
        if (h != this->marked.end()) return h;
      };
      begin = std::find_if(end, this->marked.end(), designated);
    };
    throw "Abduction error: no most plausible history amongst marked \
      histories.";
  };

  // Follows Algorithm 2 from section 3.1.2.
  // Returns true if a designated history has been added at depth
  // |sequence_of_events|.
  bool DoTL_model::abduction_step(Abduction_type abduction_type, Agent_id agent_id) {
    if (this->marked.size() == 0) {
      throw "Abduction error: no marked histories.";
    };
    std::list<History_id>::iterator marked_ptr;
    switch (abduction_type) {
      case CHRON_MIN:
        marked_ptr = this->chronological_minimization(agent_id);
        break;
      case INV_CHRON_MIN:
        marked_ptr = this->inv_chronological_minimization(agent_id);
        break;
      case UNI_EXP:
        marked_ptr = this->uniform_expansion(agent_id);
        break;
    };
    History_id h = *marked_ptr;
    this->marked.erase(marked_ptr);
    
    std::vector<History_id> current({ h });
    std::vector<History_id> children;
    std::size_t l = this->histories[h].length;
    bool added_designated_history = false;
    while (!current.empty() && l < this->sequence_of_events.size()) {
      added_designated_history = false;
      Event_model& event_model = this->sequence_of_events[l];
      
      // Histories.
      auto first_id = this->histories.size();
      for (auto& hi : current) {
        History h = this->histories[hi];
        for (auto& e: event_model.events) {
          if (this->valuate(e.precondition, h.id)) {
            auto id = this->histories.size();
            bool designated = h.designated && e.designated;
            auto he = 
              History(id, designated, h.true_propositions, h.id, h.length+1);
            he.event = e.id;
            for (const Proposition& p : e.positive_effects)
              he.true_propositions.insert(p);
            for (const Proposition& p : e.negative_effects)
              he.true_propositions.erase(p);
            add_history(he);
            added_designated_history = 
              added_designated_history || designated;
            children.emplace_back(he.id);
          };
        };
      };
      if (children.size()) l++;
      
      // Edges.
      auto last_id = this->histories.size();
      for (int agent = 0; agent < this->agents; agent++) {
        for (size_t i = first_id; i < last_id; i++) {
          this->add_plausibility_relation(agent, { {i, i} });
          if (agent == 0) this->union_relation.insert( {i, i} );
          for (auto& j : this->histories_by_length[l]) {
            if (j >= first_id && j <= i) 
              continue;
            History_id hi = this->histories[i].parent;
            History_id hj = this->histories[j].parent;
            Event_id ei = this->histories[i].event;
            Event_id ej = this->histories[j].event;
            bool hihj = 
              this->has_edge(agent, hi, hj);
            bool hjhi =
              this->has_edge(agent, hj, hi);
            bool eiej =
              event_model.has_edge(agent, ei, ej);
            bool ejei =
              event_model.has_edge(agent, ej, ei);
            if (((hihj || hjhi) && eiej && !ejei) ||
                (eiej && ejei && hihj)) {
              this->add_plausibility_relation(agent, { {i, j} });
            };
            if (((hihj || hjhi) && !eiej && ejei) ||
                (eiej && ejei && hjhi)) {
              this->add_plausibility_relation(agent, { {j, i} });
            };
            // compute union of relations.
            if (agent == 0) {
              hihj = this->union_relation.contains({hi, hj});
              hjhi = this->union_relation.contains({hj, hi});
              eiej = event_model.union_relation.contains({ei, ej});
              ejei = event_model.union_relation.contains({ej, ei});
              if (((hihj || hjhi) && eiej && !ejei) ||
                  (eiej && ejei && hihj)) {
                this->union_relation.insert({i, j});
              };
              if (((hihj || hjhi) && !eiej && ejei) ||
                  (eiej && ejei && hjhi)) {
                this->union_relation.insert({j, i});
              };
            };
          };
        }
      };

      current.clear();
      // min_children := min(i) restricted to children.
      std::vector<History_id> min_children;
      if (children.size()) {
        int num_designated = 0;
        History_id hd;
        for (auto& h : this->histories_by_length[l]) {
          if (this->is_designated(h)) {
            num_designated++;
            hd = h;
          };
        };
        if (num_designated) {
          auto hdi = this->get_information_cell(agent_id, hd);
          // Hd is closed under ~i so all histories in hdi are designated.
          if (hdi.size() == num_designated)
          {
            for (auto& h : children) {
              bool is_minimal = true;
              bool in_hdi = false;
              for (int i = 0; i < hdi.size() && is_minimal; i++) {
                in_hdi = in_hdi || hdi[i] == h;
                is_minimal = is_minimal && 
                  (hdi[i] < first_id || this->has_edge(agent_id, hdi[i], h));
              };
              if (is_minimal && in_hdi) min_children.emplace_back(h);
            };
          };
        };
      };
      
      // min_children := close(min_children) restricted to children.
      auto& tc = this->union_relation;
      bool mark;
      for (auto& h1 : children) {
        mark = true;
        for (auto& h2 : min_children) {
          if (tc.contains({h2,h1})) {
            current.emplace_back(h1);
            mark = false;
            break;
          };
        };
        if (l < this->sequence_of_events.size() && mark) {
          this->marked.push_back(h1);
        };
      };
      children.clear();
    };
    return (l == this->sequence_of_events.size()) && added_designated_history;
  };

  const std::size_t DoTL_model::size() const 
    { return this->histories.size(); };

  void DoTL_model::print() {
    for (int i = 0; i < histories_by_length.size(); i++) {
      std::cout << "Histories at level " << i << ": ";
      for (const History_id hi : this->histories_by_length[i]) {
        std::cout << "h" << this->histories[hi].id << " "; 
      };
      std::cout << std::endl;
    };
    for (const History& history : histories) {
      std::cout << "h" << history.id << " = (w" << this->histories[0].id;
      auto next = history;
      while (next.length > 0) {
        std::cout << ";" << "e" << next.event;
        next = this->histories[next.parent];
      };
      std::cout << ") :";
      for (const Proposition& p: history.true_propositions)
        std::cout << " " << p;
      std::cout << std::endl;
      for (auto agent = 0; agent < this->agents; agent++) {
        std::cout << "<=_" << agent << ": {";
        bool separator = false;
        for (const History_id most_plausible : 
            dotl::edges(this->plausibility_relation, agent, history.id)) {
          if (separator) { 
            std::cout << ",";
          }
          else {
            separator = true;
          };
          std::cout << "(h" << history.id << "," << "h" << most_plausible << ")";
        }
        std::cout << "}" << std::endl;
      }
    }
  };

  void DoTL_model::to_dot(std::unordered_map<Id, std::string> event_names,
                          std::string filename) {
    std::string colors[] =
      { "blue", "red", "green", "brown", "orange", "purple" };
    auto num_colors = sizeof(colors)/sizeof(*colors);
    std::size_t color_idx = 0;
    std::ofstream out(filename);
    out << "digraph g{" 
        << std::endl 
        << "newrank=true;"
        << std::endl 
        << "node [shape=circle, height=0.1, width=0.2, margin=\"0.5,0.5\", fixedsize=true, penwidth=0.2, fontsize=5];" 
        << std::endl 
        << "edge [fontsize=5, penwidth=0.2, arrowsize=0.3, arrowhead=open];" 
        << std::endl;

    for (auto h : this->histories) {
      std::string options = "";
      if (h.designated)
        options = "peripheries=2, ";
      //options += "xlabel=\"h";
      //options += std::to_string(h.id);
      //options += "\", ";
      options += "label=\"";
      for (auto p : h.true_propositions) {
        options += p;
        options += ",";
      }
      if (!h.true_propositions.empty())
        options.pop_back();
      options += "\"";
      out << h.id << " [" << options << "];" << std::endl;
    }
    for (auto hid : this->marked) {
      out << hid << " [style=filled, fillcolor=red];" << std::endl;
    };

    for (auto same_level_histories : histories_by_length) {
      if (!same_level_histories.empty()) {
        out << "{ rank=same; ";
        for (auto h : same_level_histories) {
          out << h << "; ";
        };
        out << "}; " << std::endl;
      };
    };

    for (auto h : this->histories) {
      bool root = h.parent == h.id;
      if (!root) {
        out << h.parent << "->" <<  h.id 
            << " [ style=dashed label=\"" << event_names[h.event]// << "(e" << h.event << ")"
            << "\" ];" << std::endl;
      };
    };

    std::unordered_map<
      std::pair<History_id, History_id>, 
      std::vector<Agent_id>,
      boost::hash<std::pair<History_id, History_id>>> to_draw;
    for (auto h : this->histories) {
      for (auto agent = 0; agent < this->agents; agent++) {
        for (const History_id most_plausible :
            this->plausibility_relation[ {agent, h.id} ].second) {
          if (h.id != most_plausible)
            to_draw[ {h.id, most_plausible} ].emplace_back(agent);
        };
      };
    };

    History_id from, to;
    std::string color;
    for (auto e : to_draw) {
      from = e.first.first;
      to = e.first.second;
      out << from << "->" << to << " [ xlabel=\"";
      for (auto i : e.second) {
        out << i << ",";
      };
      color = colors[color_idx++ % num_colors];
      out.seekp(-1, std::ios_base::end);
      out << "\", color=" <<  color << ", fontcolor=" << color << " ];" <<
        std::endl;
    }

    out << "}" << std::endl;
  };
}
