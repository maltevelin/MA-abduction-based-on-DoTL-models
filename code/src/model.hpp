#pragma once 

#include <cstddef>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <boost/variant.hpp>
#include <boost/variant/multivisitors.hpp>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <memory>
#include <random>

namespace dotl {
  
  using Id = std::size_t;
  using History_id = Id;
  using Event_id = Id;
  using Agent_id = Id;
  using Explanation = std::vector<Event_id>;
  using Proposition = char;
  using Propositions = std::unordered_set<Proposition>;

  bool in(const Propositions& P, const Proposition& p);
  
  struct Top;
  struct Bot;
  struct Neg;
  struct Con;
  struct Dis;
  struct Believes;
  struct Knows;
  struct Common_knowledge;
  struct After_event;


  struct Seq;
  struct Or;
  struct Inv;

  using Path =
    boost::variant<
      Event_id,
      boost::recursive_wrapper<Seq>,
      boost::recursive_wrapper<Or>,
      boost::recursive_wrapper<Inv>
    >;

  using Formula = 
    boost::variant<
      Proposition,
      Top,
      Bot,
      boost::recursive_wrapper<Neg>,
      boost::recursive_wrapper<Con>,
      boost::recursive_wrapper<Dis>,
      boost::recursive_wrapper<Believes>,
      boost::recursive_wrapper<Knows>,
      boost::recursive_wrapper<Common_knowledge>,
      boost::recursive_wrapper<After_event>
    >;

   struct Seq {
    Seq(const Path& pi1, const Path& pi2) : 
      pi1(std::make_shared<Path>(pi1)), pi2(std::make_shared<Path>(pi2)) {}
    const std::shared_ptr<Path> pi1, pi2;
  };
  struct Or {
    Or(const Path& pi1, const Path& pi2) : 
      pi1(std::make_shared<Path>(pi1)), pi2(std::make_shared<Path>(pi2)) {}
    const std::shared_ptr<Path> pi1, pi2;
  };
  struct Inv {
    Inv(const Path& pi) : pi(std::make_shared<Path>(pi)) {}
    const std::shared_ptr<Path> pi;
  };

  struct Top {};
  struct Bot {};
  struct Neg {
    Neg(const Formula& phi) : phi(std::make_shared<Formula>(phi)) {}
    const std::shared_ptr<Formula> phi;
  };
  struct Con {
    Con(const Formula& phi, const Formula& psi) : 
      phi(std::make_shared<Formula>(phi)), psi(std::make_shared<Formula>(psi))
    {}
    const std::shared_ptr<Formula> phi, psi;
  };
  struct Dis {
    Dis(const Formula& phi, const Formula& psi) : 
      phi(std::make_shared<Formula>(phi)), psi(std::make_shared<Formula>(psi))
    {}
    const std::shared_ptr<Formula> phi, psi;
  };
  struct Believes {
    Believes(const Formula& phi, const Agent_id& i) : 
      phi(std::make_shared<Formula>(phi)), i(i) {}
    const std::shared_ptr<Formula> phi;
    const Agent_id& i;
  };
  struct Knows {
    Knows(const Formula& phi, const Agent_id& i) : 
      phi(std::make_shared<Formula>(phi)), i(i) {}
    const std::shared_ptr<Formula> phi;
    const Agent_id& i;
  };
  struct Common_knowledge {
    Common_knowledge(const Formula& phi) : phi(std::make_shared<Formula>(phi))
    {}
    const std::shared_ptr<Formula> phi;
  };
  struct After_event {
    After_event(const Formula& phi, const Path& pi) : 
      phi(std::make_shared<Formula>(phi)), pi(std::make_shared<Path>(pi)) {}
    const std::shared_ptr<Path> pi;
    const std::shared_ptr<Formula> phi;
  };

  class DoTL_model;
  struct Valuate : public boost::static_visitor<bool> {
    bool operator() (const Proposition& p) const;
    bool operator() (const Top&) const;
    bool operator() (const Bot&) const;
    bool operator() (const Neg& phi) const;
    bool operator() (const Con& phi) const;
    bool operator() (const Dis& phi) const;
    bool operator() (const Believes& phi) const;
    bool operator() (const Knows& phi) const;
    bool operator() (const Common_knowledge& phi) const;
    bool operator() (const After_event& phi) const;
    bool operator() (const Path& pi) const;
    Valuate(const Valuate& visitor) : D(visitor.D), h(visitor.h) {}
    Valuate(DoTL_model& D, History_id h) : D(D), h(h) {}
    DoTL_model& D;
    History_id h;
  };

  struct With_path : public boost::static_visitor<bool> {
    bool operator() (const Event_id& e) const;
    bool operator() (const Seq& pi) const;
    bool operator() (const Or& pi) const;
    bool operator() (const Inv& pi) const;
    With_path(const With_path& visitor) : 
      D(visitor.D), h(visitor.h), phi(visitor.phi) {}
    With_path(const Valuate& visitor, const Formula& phi) :
      D(visitor.D), h(visitor.h), phi(phi) {}
    With_path(DoTL_model& D, History_id h, const Formula& phi) :
      D(D), h(h), phi(phi) {}
    DoTL_model& D;
    History_id h;
    const Formula& phi;
  };

  struct With_inv_path : public boost::static_visitor<bool> {
    bool operator() (const Event_id& e) const;
    bool operator() (const Seq& pi) const;
    bool operator() (const Or& pi) const;
    bool operator() (const Inv& pi) const;
    With_inv_path(const Valuate& visitor, const Formula& phi) :
      D(visitor.D), h(visitor.h), phi(phi) {}
    With_inv_path(DoTL_model& D, History_id h, const Formula& phi)
      : D(D), h(h), phi(phi) {}
    DoTL_model& D;
    History_id h;
    const Formula& phi;
  };

  struct Path_to_string : public boost::static_visitor<std::string> {
    std::string operator() (const Event_id& e) const;
    std::string operator() (const Seq& pi) const;
    std::string operator() (const Or& pi) const;
    std::string operator() (const Inv& pi) const;
  };
  std::ostream& operator<<(std::ostream&, const Path&);

  struct Formula_to_string : public boost::static_visitor<std::string> {
    std::string operator() (const Proposition& p) const;
    std::string operator() (const Top&) const;
    std::string operator() (const Bot&) const;
    std::string operator() (const Neg&) const;
    std::string operator() (const Con& phi) const;
    std::string operator() (const Dis& phi) const;
    std::string operator() (const Believes& phi) const;
    std::string operator() (const Knows& phi) const;
    std::string operator() (const Common_knowledge& phi) const;
    std::string operator() (const After_event& phi) const;
  };
  std::ostream& operator<<(std::ostream&, const Formula&);

  class World {
    public:
      World(Id id) :
        id(id), designated(false) {}
      World(Id id, bool designated) :
        id(id), designated(designated) {}

      Id id;
      bool designated;
  };

  class Event : public World {
    public:
      Event(Event_id id, bool designated,
            Formula precondition,
            std::initializer_list<Proposition> positive_effects,
            std::initializer_list<Proposition> negative_effects) :
        World(id, designated),
        positive_effects(Propositions(positive_effects)),
        negative_effects(Propositions(negative_effects)),
        precondition(precondition) {}

      Event copy(const Event& event) const;
      bool operator!=(const Event& other) const;

      Propositions positive_effects;
      Propositions negative_effects;
      Formula precondition;
  };

  using Edge = std::pair<Id, Id>;
  bool operator<(const Edge& left, const Edge& right);
  using Relation = 
    std::unordered_map<std::pair<Agent_id, Id>, 
                       std::pair<bool, std::vector<Id>>,
                       boost::hash<std::pair<size_t, size_t>>>;

  std::vector<Id>& edges(Relation& R, Agent_id i, Id from);

  const bool has_edge(Relation& R, Agent_id i, Id from, Id to);

  void add_edges_to_relation(Relation& relation, Agent_id agent_id,
      std::vector<Edge> edges);
  void add_inverse_edges_to_relation(Relation& relation, Agent_id agent_id,
      std::vector<Edge> edges);

  class Event_model {
    public:
      Event_model(std::size_t agents) : agents(agents), offset(total_events) {}

      const Propositions& get_positive_effects(Event_id event_id) const;
      const Propositions& get_negative_effects(Event_id event_id) const;
      const Formula& get_precondition(Event_id event_id) const;
      bool is_designated(Event_id) const;
      bool valuate(const Formula&, Event_id event_id) const;
      const bool has_edge(Agent_id i, Event_id from, Event_id to);
      const bool indistinguishable_from(Agent_id i, Event_id e1, Event_id e2);

      Event_id new_event(bool designated, 
          Formula& precondition,
          std::initializer_list<Proposition> positive_effects,
          std::initializer_list<Proposition> negative_effects);
      void add_plausibility_relation(Agent_id agent_id, std::vector<Edge>
          edges);
      void perspective_shift(Agent_id agent_id);

      void print();

      std::vector<Event> events;
      std::unordered_set<Edge, boost::hash<std::pair<size_t, size_t>>>
        union_relation;
    private:
      static Event_id total_events;
      const Event_id offset;
      Relation plausibility_relation;
      std::size_t agents;
  };

  class History : public World {
    public:
      History(History_id id) : 
        World(id), true_propositions(Propositions()), length(0) {}
      History(History_id id, bool designated, Propositions true_propositions,
          History_id parent, size_t length) :
        World(id, designated), true_propositions(true_propositions),
        parent(parent), length(length) {}

      History copy(const History& history) const;
      bool operator<(const History& other) const;

      Propositions true_propositions;
      size_t length;
      Event_id event;
      History_id parent;
  };

  enum Abduction_type { CHRON_MIN, INV_CHRON_MIN, UNI_EXP }; 

  class DoTL_model {
    public:
      // This constructor implies that there can be only one initial history.
      // Uncertainty in initial situations is added by a product update.
      DoTL_model(std::size_t agents, std::initializer_list<Proposition>
          true_propositions);

      const Propositions& get_true_propositions(History_id history_id) const;

      std::vector<History_id> get_group_reachable(History_id history_id);
      std::vector<History_id> get_information_cell(Agent_id agent_id,
          History_id history_id);
      std::vector<History_id> get_minimal_histories(Agent_id agent_id,
          History_id history_id);
      std::vector<History_id> get_minimal_histories(Agent_id agent_id,
          std::vector<History_id>& information_cell);
      std::vector<History_id>& edges_from(Agent_id agent_id, History_id
          history_id);
      std::vector<History_id>& edges_to(Agent_id agent_id, History_id
          history_id);
      const std::optional<History_id> get_child(History_id history_id, Event_id
          event_id) const;
      const std::optional<History_id> get_parent(History_id history_id, Event_id
          event_id) const;

      bool is_designated(History_id) const;
      bool valuate(const Formula& phi, History_id history_id);
      bool valuate(const Formula& phi);

      const bool has_edge(Agent_id i, History_id from, History_id to);
      const bool indistinguishable_from(
          Agent_id i, History_id e1, History_id e2);

      // restricted_update(E, NULL) = restricted_update(E, frontier(D)).
      bool restricted_update(
          Event_model& event_model, std::vector<History_id>* I);

      void action_priority_update(Event_model& event_model);
      bool most_plausible_update(Event_model& event_model, Agent_id agent_id);
      bool abduction_step(Abduction_type abduction_type, Agent_id agent_id);

      const std::size_t size() const;
      void print();
      void to_dot(std::unordered_map<Id, std::string> event_names, 
                  std::string filename);

    private:
      void add_plausibility_relation(Agent_id agent_id, std::vector<Edge>
          edges);
      std::size_t agents;
      bool find_history(std::vector<History>& hs, History& h);
      bool find_history_id(std::vector<History_id>& hs, History_id& h);
      std::vector<History_id>& frontier();
      void add_history(const History& h);
      std::size_t max_length;
      std::vector<Event_model> sequence_of_events;
      // histories[h.id] = h.
      std::vector<History> histories;
      // in general, histories[h.length][h.id] != h.
      std::vector<std::vector<History_id>> histories_by_length;
      std::list<History_id> marked;
      Relation plausibility_relation;
      Relation inverse_plausibility_relation;
      std::unordered_set<Edge, boost::hash<std::pair<size_t, size_t>>>
        union_relation;
      std::list<History_id>::iterator most_plaus_marked_from_sorted(
          Agent_id agent_id);
      std::list<History_id>::iterator chronological_minimization(
          Agent_id agent_id);
      std::list<History_id>::iterator inv_chronological_minimization(
          Agent_id agent_id);
      std::list<History_id>::iterator uniform_expansion(
          Agent_id agent_id);
      std::default_random_engine rng;
  };
}
